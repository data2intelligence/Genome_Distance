from django.contrib import admin

from .models import AliasName, Gene, Probe, GeneProbe, Loci

class AliasNameAdmin(admin.ModelAdmin):
    ordering = ('name',)
    list_display = ['name', 'taxid']
    search_fields = ('name', 'taxid')

admin.site.register(AliasName, AliasNameAdmin)


class GeneAdmin(admin.ModelAdmin):
    ordering = ('Symbol',)
    list_display = ['Symbol', 'taxid', 'get_alias']
    search_fields = ('Symbol', 'taxid', 'Alias__name')
    
    def get_alias(self, obj):
        return ' | '.join([p.name for p in obj.Alias.all()])

admin.site.register(Gene, GeneAdmin)

    
class ProbeAdmin(admin.ModelAdmin):
    ordering = ('ID',)
    list_display = ['ID', 'Chrom', 'Start', 'End', 'taxid']
    search_fields = ('ID', 'Chrom', 'Start', 'End', 'taxid')

admin.site.register(Probe, ProbeAdmin)


class GeneProbeAdmin(admin.ModelAdmin):
    list_display = ('Gene', 'Probe', 'distance')
    ordering = ('distance',)
    search_fields = ('Gene__Symbol', 'Probe__ID')

admin.site.register(GeneProbe, GeneProbeAdmin)


class LociAdmin(admin.ModelAdmin):
    list_display = ('Probe', 'Cell', 'Condition', 'x', 'y', 'z')
    ordering = ('Probe','Cell', 'Condition', )
    search_fields = ('Probe__ID', 'Cell', 'Condition')

admin.site.register(Loci, LociAdmin)
