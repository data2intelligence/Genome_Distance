from django.db import models
from Genome_Distance.settings import MAX_CHAR_LENGTH

# Create your models here.
class AliasName(models.Model):
    name = models.CharField(max_length=MAX_CHAR_LENGTH, primary_key = True)
    taxid = models.IntegerField()
    
    def __unicode__(self): return self.name
    def __str__(self): return self.__unicode__()


class Gene(models.Model):
    Symbol = models.CharField(max_length=MAX_CHAR_LENGTH, primary_key = True)
    taxid = models.IntegerField()
    
    # transcription start site
    Chrom = models.CharField(max_length=MAX_CHAR_LENGTH, null=True)
    TSS = models.IntegerField(null=True)
    
    # use for searching
    Alias = models.ManyToManyField(AliasName)
    
    def __unicode__(self): return self.Symbol
    
    def __str__(self): return self.__unicode__()



class Probe(models.Model):
    ID = models.CharField(max_length=MAX_CHAR_LENGTH, primary_key = True)
    Chrom = models.CharField(max_length=MAX_CHAR_LENGTH)
    Start = models.IntegerField()
    End = models.IntegerField()
    taxid = models.IntegerField()
    
    def __unicode__(self): return '%s %s (%d, %d) %d' % (self.ID, self.Chrom, self.Start, self.End, self.taxid)
    def __str__(self): return self.__unicode__()


class GeneProbe(models.Model):
    Gene = models.ForeignKey(Gene, related_name='geneprobe', on_delete=models.CASCADE)
    Probe = models.ForeignKey(Probe, related_name='geneprobe', on_delete=models.CASCADE)
    
    distance = models.IntegerField()
    
    class Meta:
        unique_together = (('Gene', 'Probe'),)
    
    def __unicode__(self): return '%s %s %d' % (self.Gene.Symbol, self.Probe.ID, self.distance)
    def __str__(self): return self.__unicode__()


class Loci(models.Model):
    Probe = models.ForeignKey(Probe, related_name='loci', on_delete=models.CASCADE)
    
    Cell = models.CharField(max_length=MAX_CHAR_LENGTH)
    Condition = models.CharField(max_length=MAX_CHAR_LENGTH)
    
    x = models.FloatField()
    y = models.FloatField()
    z = models.FloatField()
    
    class Meta:
        unique_together = (('Probe', 'Cell', 'Condition'),)
    
    def __unicode__(self): return '%s %s %s (%.2f, %.2f, %.2f)' % (self.Probe.ID, self.Cell, self.Condition, self.x, self.y, self.z)
    def __str__(self): return self.__unicode__()


class Background(models.Model):
    Condition = models.CharField(max_length=MAX_CHAR_LENGTH)
    
    Chrom_A = models.CharField(max_length=MAX_CHAR_LENGTH)
    Chrom_B = models.CharField(max_length=MAX_CHAR_LENGTH)
    
    histogram_mean = models.JSONField()
    histogram_std = models.JSONField()
    
    class Meta:
        unique_together = (('Chrom_A', 'Chrom_B', 'Condition'),)
    
    def __unicode__(self):
        return '%s %s %s %s %s' % (self.Condition, self.Chrom_A, self.Chrom_B, self.histogram_mean, self.histogram_std)
    
    def __str__(self): return self.__unicode__()
