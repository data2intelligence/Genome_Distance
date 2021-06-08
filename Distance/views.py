from django.shortcuts import render
from Distance.models import Gene, AliasName, Probe, GeneProbe, Loci, Background
from Genome_Distance.settings import HALF_DISTANCE, HISTO_START, HISTO_END, HISTO_STEP

import pandas, numpy

histo_bins = numpy.arange(HISTO_START, HISTO_END + HISTO_STEP, HISTO_STEP)


def search_gene(gene):
    gene_lst = Gene.objects.filter(Symbol=gene)
    if len(gene_lst) > 0: return gene_lst
    
    alias_lst = AliasName.objects.filter(name=gene)
    if len(alias_lst) > 0:
        genes = set()
        for alias in alias_lst:
            genes.update(alias.gene_set.all())
        return genes
    
    return None


def search_probe(gene, half_distance):
    # case 1: gene 5' probe is already available
    probe_lst = Probe.objects.filter(ID=gene.Symbol)
    
    if probe_lst.count() == 1:
        return [[probe_lst[0], 1]]
    else:
        assert probe_lst.count() == 0
    
    # case 2: find nearest probe to 5'
    probe_lst = GeneProbe.objects.filter(Gene=gene)
    
    a = numpy.log(2)/half_distance
    
    if probe_lst.count() > 0:
        results = []
        
        for probe in probe_lst:
            results.append([probe.Probe, numpy.exp(-a*probe.distance)])
        
        return results
    
    return None



def get_distance(loci_lst_a, loci_lst_b):
    common = set([loci.Cell + '\t' + loci.Condition for loci in loci_lst_a])
    common.intersection_update([loci.Cell + '\t' + loci.Condition for loci in loci_lst_b])
    
    condition_map = {}
    
    for entry in common:
        cell, condition = entry.split('\t',1)
        
        loci_a = loci_lst_a.filter(Cell=cell, Condition=condition)
        assert loci_a.count() == 1
        loci_a = loci_a[0]
        
        loci_b = loci_lst_b.filter(Cell=cell, Condition=condition)
        assert loci_b.count() == 1
        loci_b = loci_b[0]
        
        dis = (103*(loci_a.x - loci_b.x))**2 + (103*(loci_a.y - loci_b.y))**2 + (250*(loci_a.z - loci_b.z))**2
        dis = numpy.sqrt(dis)/1000  # nM -> uM
        
        lst = condition_map.get(condition)
        if lst is None:
            lst = condition_map[condition] = []
        lst.append([cell, dis])
    
    return condition_map
        


def json_to_series(arr, title):
    arr = pandas.read_json(arr, typ='series', convert_axes=False)
    arr.index = arr.index.astype(float)
    arr.name = title
    
    return arr     



def search(request):
    select_a = select_b = query_a = query_b = query_a_lst = query_b_lst = None
    
    if 'query_input_a' in request.POST:
        query_a = request.POST['query_input_a']
        
        if not query_a:
            return render(request, 'error.html', {'message': 'Please enter a search term for gene A.'})
        
        query_a_lst = search_gene(query_a)
        
        if not query_a_lst:
            return render(request, 'error.html', {'message': 'Cannot find gene for query A %s.' % query_a})

    if 'query_input_b' in request.POST:
        query_b = request.POST['query_input_b']
        
        if not query_b:
            return render(request, 'error.html', {'message': 'Please enter a search term for gene B.'})
        
        query_b_lst = search_gene(query_b)
        
        if not query_b_lst:
            return render(request, 'error.html', {'message': 'Cannot find gene for query B %s.' % query_b})
    
    
    if query_a_lst and query_b_lst:
        
        return render(request, 'search_results.html', {
            'query_a': query_a,
            'query_a_lst': query_a_lst,
            
            'query_b': query_b,
            'query_b_lst': query_b_lst,
            
            'half_distance': HALF_DISTANCE,
            })
    
    
    
    if 'select_a' in request.GET:
        select_a = request.GET['select_a']
        
    if 'select_b' in request.GET:
        select_b = request.GET['select_b']
        
    
    if select_a and select_b:
        gene_a = Gene.objects.get(Symbol=select_a)
        gene_b = Gene.objects.get(Symbol=select_b)
        
        select_a = gene_a
        select_b = gene_b
        
        if gene_a.taxid != gene_b.taxid:
            return render(request, 'error.html', {
                'message': 'Different Tax IDs %d, %d between %s and %s.' % (gene_a.taxid, gene_b.taxid, gene_a.Symbol, gene_b.Symbol)
                })
        
        half_distance = float(request.GET['half_distance'])
        
        probe_lst_a = search_probe(gene_a, half_distance)
        if probe_lst_a is None: return render(request, 'error.html', {'message': 'Cannot find probes for %s.' % gene_a.Symbol})
        
        probe_lst_b = search_probe(gene_b, half_distance)
        if probe_lst_b is None: return render(request, 'error.html', {'message': 'Cannot find probes for %s.' % gene_b.Symbol})
        
        result_lst = []
        for probe_a, confi_a in probe_lst_a:
            loci_lst_a = Loci.objects.filter(Probe=probe_a)
            
            for probe_b, confi_b in probe_lst_b:
                loci_lst_b = Loci.objects.filter(Probe=probe_b)
                
                condition_map = get_distance(loci_lst_a, loci_lst_b)
                
                for condition, distance_table in condition_map.items():
                    background = Background.objects.filter(Condition=condition, Chrom_A=probe_a.Chrom, Chrom_B=probe_b.Chrom)
                    
                    if background.count() == 0:
                        background = Background.objects.filter(Condition=condition, Chrom_A=probe_b.Chrom, Chrom_B=probe_a.Chrom)
                    
                    background_bins = background_mean = background_std = None
                    
                    if background.count() > 0:
                        assert background.count() == 1
                        
                        background_mean = json_to_series(background[0].histogram_mean, 'Mean')
                        background_std = json_to_series(background[0].histogram_std, 'Std')
                        background_bins = background_mean.index
                        
                        assert sum(background_bins != background_std.index) == 0
                        
                        # modal normalization
                        background_mean /= numpy.max(background_mean)
                    
                    # generate histogram normalized to fractions
                    histo = numpy.histogram([min(v[1], HISTO_END) for v in distance_table], bins=histo_bins)[0] / len(distance_table)
                    
                    # modal normalization
                    histo /= numpy.max(histo)
                    
                    result_lst.append([probe_a, confi_a, probe_b, confi_b, condition, distance_table, 
                                       histo_bins[:-1], histo,
                                       background_bins, background_mean])
        
        return render(request, 'distance.html', {
                'select_a': select_a,
                'select_b': select_b,
                'result_lst': result_lst,
                'bar_width': HISTO_STEP,
                })
    
    # {'genes': [v.Symbol for v in Gene.objects.all()]}
    return render(request, 'search.html')
