import os, sys, pandas
from django.core.management.base import BaseCommand#, CommandError
from Distance.models import AliasName, Gene, Probe, GeneProbe, Loci


class Command(BaseCommand):
    help = 'Build the database data from command line'
    
    def add_arguments(self, parser):
        parser.add_argument('operation', nargs='+', type=str)
    
    def handle(self, *args, **options):
        operation = options['operation'][0]
        
        if operation == 'insert_gene':
            taxid = int(options['operation'][1])
            self.insert_gene(taxid)
        
        elif operation == 'insert_probe':
            coordinate = options['operation'][1]
            taxid = int(options['operation'][2])
            self.insert_probe(coordinate, taxid)
        
        elif operation == 'insert_geneprobe':
            self.insert_geneprobe(options['operation'][1])
        
        elif operation == 'insert_loci':
            condition = options['operation'][1]
            data = options['operation'][2]
            
            self.insert_loci(condition, data)
        
        else:
            sys.stderr.write('Cannot recognize operation %s\n' % operation)
        
    
    def insert_gene(self, taxid):
        cnt_genes = cnt_alias = cnt_previous = 0
        
        gene_info = os.path.join(os.getenv("HOME"), 'workspace', 'Data', 'GeneInfo', 'NCBI', 'gene_info.%d' % taxid)
        
        fin = open(gene_info)
        
        for l in fin:
            fields = l.rstrip('\n').split('\t')
            
            symbol = fields[2].strip()
            alias = fields[4].strip()
            
            if len(symbol) == 0: continue
            
            rc_gene = Gene.objects.filter(Symbol=symbol)
            
            if rc_gene.count() == 0:
                rc_gene = Gene.objects.create(Symbol=symbol, taxid=taxid)
                cnt_genes += 1
            else:
                assert rc_gene.count() == 1
                rc_gene = rc_gene[0]
                cnt_previous += 1
            
            # insert all alias if there is any
            if len(alias) == 0: continue
            
            alias_lst = alias.split('|')
            
            for alias in alias_lst:
                if len(alias) == 0: continue
                
                rc_alias = AliasName.objects.filter(name=alias)
                
                if rc_alias.count() == 0:
                    rc_alias = AliasName.objects.create(name=alias, taxid=taxid)
                    cnt_alias = cnt_alias + 1
                else:
                    assert rc_alias.count() == 1
                    rc_alias = rc_alias[0]
                
                rc_gene.Alias.add(rc_alias)
        
        fin.close()
        
        print(cnt_genes, 'created genes', cnt_previous, 'previous genes', cnt_alias, 'alias')
    
    
    def insert_probe(self, coordinate, taxid):
        # Step 1: read probe coordinates
        coordinate = pandas.read_excel(coordinate, engine='openpyxl')
        assert coordinate['Name'].value_counts().max() == 1
        
        cnt_probe = cnt_previous = 0
        
        for _, arr in coordinate.iterrows():
            rc_probe = Probe.objects.filter(ID = arr['Name'])
            
            if rc_probe.count() == 0:
                rc_probe = Probe.objects.create(ID=arr['Name'], Chrom=arr['Chrom'], Start=arr['Start'], End=arr['End'], taxid=taxid)
                cnt_probe += 1
            else:
                assert rc_probe.count() == 1
                cnt_previous += 1
            
        print(cnt_probe, 'created probed', cnt_previous, 'previous probes')
        
    
    def insert_geneprobe(self, f):
        data = pandas.read_csv(f, sep='\t')
        
        cnt_geneprobe = cnt_previous = 0
        
        for _, arr in data.iterrows():
            gene = Gene.objects.filter(Symbol=arr['Gene'])
            if gene.count() == 0: continue
            
            probe = Probe.objects.filter(ID=arr['Probe'])
            if probe.count() == 0: continue
            
            assert gene.count() == 1 and probe.count() == 1
            gene = gene[0]
            probe = probe[0]
            
            rc_geneprobe = GeneProbe.objects.filter(Gene=gene, Probe=probe)
            
            if rc_geneprobe.count() == 0:
                rc_geneprobe = GeneProbe.objects.create(Gene=gene, Probe=probe, distance=arr['Distance'])
                cnt_geneprobe += 1
            else:
                assert rc_geneprobe.count() == 1
                cnt_previous += 1
                    
        print(cnt_geneprobe, 'created geneprobe', cnt_previous, 'previous geneprobe')

    

    def insert_loci(self, condition, data):
        data = pandas.read_csv(data)
    
        # cell ID is fov(field of view) plus cell ID within each fov
        flag = data.loc[:, 'fov'].astype(str) + '_' + data.loc[:, 'cellID'].astype(str) + '\t' + data.loc[:, 'geneID']
        data = data.loc[:, ['x', 'y', 'z']].groupby(flag).median()
        data[['cell', 'probe']] = [v.split('\t',1) for v in data.index]
        
        cnt_loci = cnt_previous = 0
        
        for _, arr in data.iterrows():
            probe = Probe.objects.filter(ID=arr['probe'])
                
            if probe.count() == 0:
                sys.stderr.write('Cannot find probe %s\n' % arr['probe'])
                continue
            
            assert probe.count() == 1
            probe = probe[0]
            
            rc_loci = Loci.objects.filter(Probe=probe, Cell=arr['cell'], Condition=condition)
            
            if rc_loci.count() == 0:
                rc_loci = Loci.objects.create(Probe=probe, Cell=arr['cell'], Condition=condition, x=arr['x'], y=arr['y'], z=arr['z'])
                cnt_loci += 1
            else:
                assert rc_loci.count() == 1
                cnt_previous += 1
        
        print(cnt_loci, 'created loci', cnt_previous, 'previous loci')
