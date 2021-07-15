# vim: syntax=python

# >6.0 is required for modules
from snakemake.utils import min_version
min_version("6.0")

import pandas as pd

chrs = [ x for x in range(1, 23) ]   # Chromosomes to analyze (autosomes)
output_plot_and_table = [ 'svg', 'pdf', 'tsv' ]  # Common plot+table output

# Read in the input manifest
manifest = pd.read_csv('/n/data1/hms/dbmi/park/jluquette/glia/analysis/INPUT_MANIFEST',
    sep='\t',comment='#')
celltypes = set(manifest['celltype'])
qualtypes = set(manifest['qualtype'])
celltypes_to_compute = [ 'neuron', 'oligo' ]
snv_qualtypes = [ 'A', 'AB', 'ABM' ]
indel_qualtypes = [ 'indel_A' ]
all_qualtypes = snv_qualtypes + indel_qualtypes

# Sample-specific metadata
# XXX: FIXME: these metadata files are actually produced by this snakefile.
path_to_metadata = "/n/data1/hms/dbmi/park/jluquette/glia/figures/test2/input"

neuron_manifest = pd.read_csv(path_to_metadata + "/neuron___meta___A.csv")
neuron_donors = [ str(x) for x in neuron_manifest['donor'].to_list() ]
neuron_samples = neuron_manifest['sample'].to_list()
neuron_paths = [ '/n/data1/hms/dbmi/park/jluquette/pta/' + neuron_donors[i] + '/scansnv_fdr01_noX/callable_regions/' + neuron_samples[i] for i in range(0, len(neuron_samples)) ]
neuron_dict = dict(zip(neuron_samples, neuron_paths))


oligo_manifest = pd.read_csv(path_to_metadata + "/oligo___meta___A.csv")
oligo_donors = [ str(x) for x in oligo_manifest['donor'].to_list() ]
oligo_samples = oligo_manifest['sample'].to_list()
oligo_paths = [ '/n/data1/hms/dbmi/park/jluquette/glia/' + oligo_donors[i] + '/scansnv/callable_regions/' + oligo_samples[i] for i in range(0, len(oligo_samples)) ]
oligo_dict = dict(zip(oligo_samples, oligo_paths))

pathdict = { **neuron_dict, **oligo_dict }


# Roadmap Epigenomics metadata
#roadmap_metadata = pd.read_csv(config['ROADMAP_EIDS'])
#roadmap_eids = roadmap_metadata['EID']
#roadmap_eids = [ 'E073' ]
#roadmap_marks = [ 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3' ]
binsizes = [ '1000', '10000', '100000', '1000000' ]
qsizes = [ 3, 5, 10, 50 ]


# Constant config variables for the enrichment pipeline

enrichment_config = {
    # The otput directory and signal files are the only parameters that should
    # change between analyses.
    #'output_dir': 'path/to/output',
    #'SIGNAL_MANIFEST': 'path/to/sig.csv',

    # These parameters should be constant
    'qbed_from_bigwig_script': '/n/data1/hms/dbmi/park/jluquette/glia/analysis/scripts/make_qbed_from_bigwig.sh',
    'MUT_MANIFEST': '/n/data1/hms/dbmi/park/jluquette/glia/analysis/MUTATION_MANIFEST',
    'quantiles': qsizes,
    # must match binsizes
    'tiles': {
        '1000000': 'alignability/genome_tiles/genome_tiles_1000000binsize.bed',
        '100000':  'alignability/genome_tiles/genome_tiles_100000binsize.bed',
        '10000':   'alignability/genome_tiles/genome_tiles_10000binsize.bed',
        '1000':    'alignability/genome_tiles/genome_tiles_1000binsize.bed'
    },
    'mut_to_perm': {
        'neuron___A': 'input/neuron___perm___A.rda',
        'neuron___AB': 'input/neuron___perm___AB.rda',
        'neuron___ABM': 'input/neuron___perm___ABM.rda',
        'neuron___indel_A': 'input/neuron___perm___indel_A.rda',
        'oligo___A': 'input/oligo___perm___A.rda',
        'oligo___AB': 'input/oligo___perm___AB.rda',
        'oligo___ABM': 'input/oligo___perm___ABM.rda',
        'oligo___indel_A': 'input/oligo___perm___indel_A.rda'
    }
}

wildcard_constraints:
    celltype='|'.join(celltypes),
    collapsed='|COLLAPSED.',
    qualtype='|'.join(qualtypes),
    cosmic='cosmic_full|cosmic_reduced',
    umap='atac|rna',
    datasource="roadmap|encode/replication_timing"  # still used?


rule all:
    input:
        # alignability
        expand("alignability/plots/chromosome_bin_classes_heatmap.{output}",
            output=['svg', 'pdf', 'jpeg']),
        expand("alignability/plots/chromosome_bin_classes_barplot.{output}",
            output=['svg', 'pdf']),
        expand("alignability/data_plots/chr{chr}.{resolution}.{output}",
            chr=chrs, resolution=[ '1k', '10k', '100k', '1m' ],
            output=['svg', 'pdf']),
        # 1b
        expand('plots/fig1/distribution_analysis/{qualtype}_gene_classes.{collapsed}{output}',
            qualtype=snv_qualtypes,  # indels are always added
            collapsed=[ '', 'COLLAPSED.'],
            output=output_plot_and_table),
        # 1c
        expand("plots/fig1/indel_size_analysis/indel_size.{output}",
            output=[ 'pdf', 'svg']),
        # 1d
        expand("plots/fig1/snpeff_analysis/snpeff_analysis.{output}",
            output=output_plot_and_table),
        # 1e,f
        expand('plots/fig1/aging_rate_analysis/{qualtype}.{output}',
            qualtype=[ 'A', 'indel_A' ],
            output=[ 'pdf', 'svg']),
        expand('plots/fig1/aging_rate_analysis/{qualtype}_{table}.csv',
            qualtype=[ 'A', 'indel_A' ],
            table=[ 'burdens', 'model' ]),
        # 2a
        expand("plots/fig2/raw_spectrum/snv_spectrum.{output}",
            output=output_plot_and_table),
        # 2b
        expand("analysis/fig2/cosmic_signature_inclusion/{celltype}_signature_scores_{qualtype}.csv",
            celltype=celltypes_to_compute,
            qualtype=[ 'A', 'indel_A' ]),
        # probably supplementary, related to fig2
        expand("analysis/fig2/cosmic_signature_inclusion/final_signature_selection_{qualtype}.{output}",
            qualtype=[ 'A', 'indel_A' ],
            output=[ 'pdf', 'svg', 'csv' ]),
        # 2b,f
        expand("plots/fig2/cosmic_aging/{celltype}_barplots_{qualtype}.{output}",
            celltype=celltypes_to_compute,
            qualtype=[ 'A', 'indel_A' ],
            output=[ 'pdf', 'svg' ]),
        # 2c,g
        expand("plots/fig2/cosmic_aging/scatterplots_{qualtype}.{output}",
            qualtype=[ 'A', 'indel_A' ],
            output=[ 'pdf', 'svg', 'csv' ]),
        # 2e
        expand("plots/fig2/raw_spectrum/indel_spectrum.{output}",
            output=output_plot_and_table),
        # 3a
        expand("plots/fig3/scrnaseq/umap_plot.{output}",
            output=[ 'svg', 'pdf', 'jpeg' ]),
        # 3b
        expand("plots/fig3/scatacseq/umap_plot.{output}",
            output=[ 'svg', 'pdf', 'jpeg' ]),
        # 3c
        expand('plots/fig3/encode_replichip/quantile/{celltype}___{qualtype}.{nquantiles}quantiles.{output}',
            celltype=celltypes_to_compute,
            qualtype=all_qualtypes,
            nquantiles=qsizes,
            output=[ 'svg', 'pdf' ]),
        # 3e
        expand('plots/fig3/roadmap_enrichment/quantile/{celltype}___{qualtype}.{nquantiles}quantiles.{output}',
            celltype=celltypes_to_compute,
            qualtype=all_qualtypes,
            nquantiles=qsizes,
            output=[ 'svg', 'pdf' ]),
        expand('plots/enrichment/{datasource}/quantile/{celltype}___{qualtype}.{nquantiles}quantiles.{output}',
            celltype=celltypes_to_compute,
            qualtype=all_qualtypes,
            nquantiles=qsizes,
            datasource=[ 'conservation', 'boca', 'depth', 'nott', 'rizzardi' ],
            #datasource=[ 'rizzardi' ],
            output=[ 'svg', 'pdf', 'csv' ])
        # Just signals for now, not sure which figure (if any) will use them
        #expand('enrichment/{datasource}/quantile/{celltype}___{qualtype}/{nquantiles}quantiles.csv',
            #celltype=celltypes_to_compute,
            #qualtype=all_qualtypes,
            #nquantiles=qsizes,
            #datasource=[ 'conservation', 'boca', 'depth', 'nott', 'rizzardi' ])


include: "snakefile.data"
include: "snakefile.alignability"
include: "snakefile.fig1"
include: "snakefile.fig2"
include: "snakefile.fig3"


# Below is several instances of
# boilerplate code to use the enrichment module to automatically
# run bigWig signal files through the enrichment pipeline.

########################################################################
# Roadmap epigenomics histones
########################################################################
enrichment_roadmap_config = dict(
    **{ 'output_dir': 'enrichment/roadmap',
        'SIGNAL_MANIFEST': '/n/data1/hms/dbmi/park/jluquette/glia/analysis/ROADMAP_HISTONE.MANIFEST' },
    **enrichment_config
)

module enrichment_roadmap:
    snakefile: "snakefile.enrichment"
    config: enrichment_roadmap_config

use rule * from enrichment_roadmap as enrichment_roadmap_*

use rule enrichment_plot from enrichment_roadmap as enrichment_roadmap_enrichment_plot with:
    params:
        ignore='enrichment_grid_R_ignore_none',
        group='mark'
    output:
        expand('plots/enrichment/roadmap/quantile/{{mutclass}}.{{nquantiles}}quantiles.{output}',
            output=[ 'svg', 'pdf', 'csv' ])
    log:
        'plots/enrichment/roadmap/quantile/{mutclass}.{nquantiles}quantiles.log'


########################################################################
# ENCODE replication timing as measured by Repli-chip
########################################################################
enrichment_replichip_config = dict(
    **{ 'output_dir': 'enrichment/encode_replichip',
        'SIGNAL_MANIFEST': '/n/data1/hms/dbmi/park/jluquette/glia/analysis/ENCODE_REPLICHIP.MANIFEST' },
    **enrichment_config
)

module enrichment_replichip:
    snakefile: "snakefile.enrichment"
    config: enrichment_replichip_config

use rule * from enrichment_replichip as enrichment_replichip_*

use rule enrichment_plot from enrichment_replichip as enrichment_replichip_enrichment_plot with:
    params:
        ignore='BINSIZE=1000',
        group='datasource'
    output:
        expand('plots/enrichment/replichip/quantile/{{mutclass}}.{{nquantiles}}quantiles.{output}',
            output=[ 'svg', 'pdf', 'csv' ])
    log:
        'plots/enrichment/replichip/quantile/{mutclass}.{nquantiles}quantiles.log'


########################################################################
# Conservation tracks from UCSC
########################################################################
enrichment_conservation_config = dict(
    **{ 'output_dir': 'enrichment/conservation',
        'SIGNAL_MANIFEST': '/n/data1/hms/dbmi/park/jluquette/glia/analysis/CONSERVATION.MANIFEST' },
    **enrichment_config
)

module enrichment_conservation:
    snakefile: "snakefile.enrichment"
    config: enrichment_conservation_config

use rule * from enrichment_conservation as enrichment_conservation_*

# Conservation tracks have bp resolution over the entire genome. This
# causes bigWigAverageOverBed to consume unbelievably large amounts of
# memory for large windows like 1MB.
use rule make_qbed_from_bigwig from enrichment_conservation as enrichment_conservation_make_qbed_from_bigwig with:
    resources:
        mem=65000

use rule enrichment_plot from enrichment_conservation as enrichment_conservation_enrichment_plot with:
    params:
        ignore='enrichment_grid_R_ignore_none',
        group='track'
    output:
        expand('plots/enrichment/conservation/quantile/{{mutclass}}.{{nquantiles}}quantiles.{output}',
            output=[ 'svg', 'pdf', 'csv' ])
    log:
        'plots/enrichment/conservation/quantile/{mutclass}.{nquantiles}quantiles.log'


########################################################################
# BOCA brain region ATAC-seq
########################################################################
enrichment_boca_config = dict(
    **{ 'output_dir': 'enrichment/boca',
        'SIGNAL_MANIFEST': '/n/data1/hms/dbmi/park/jluquette/glia/analysis/BOCA.MANIFEST' },
    **enrichment_config
)

module enrichment_boca:
    snakefile: "snakefile.enrichment"
    config: enrichment_boca_config

use rule * from enrichment_boca as enrichment_boca_*

# BOCA ATAC tracks are ~10-fold larger than most others, so give them a
# bit more RAM.
use rule make_qbed_from_bigwig from enrichment_boca as enrichment_boca_make_qbed_from_bigwig with:
    resources:
        mem=15000

use rule enrichment_plot from enrichment_boca as enrichment_boca_enrichment_plot with:
    params:
        ignore='enrichment_grid_R_ignore_none',
        group='datasource,celltype,region'
    output:
        expand('plots/enrichment/boca/quantile/{{mutclass}}.{{nquantiles}}quantiles.{output}',
            output=[ 'svg', 'pdf', 'csv' ])
    log:
        'plots/enrichment/boca/quantile/{mutclass}.{nquantiles}quantiles.log'


########################################################################
# Sequncing depth
########################################################################
enrichment_depth_config = dict(
    **{ 'output_dir': 'enrichment/depth',
        'SIGNAL_MANIFEST': '/n/data1/hms/dbmi/park/jluquette/glia/analysis/DEPTH.MANIFEST' },
    **enrichment_config
)

module enrichment_depth:
    snakefile: "snakefile.enrichment"
    config: enrichment_depth_config

use rule * from enrichment_depth  as enrichment_depth_*

use rule enrichment_plot from enrichment_depth as enrichment_depth_enrichment_plot with:
    params:
        ignore='enrichment_grid_R_ignore_none',
        group='celltype'
    output:
        expand('plots/enrichment/depth/quantile/{{mutclass}}.{{nquantiles}}quantiles.{output}',
            output=[ 'svg', 'pdf', 'csv' ])
    log:
        'plots/enrichment/depth/quantile/{mutclass}.{nquantiles}quantiles.log'


########################################################################
# Nott et al histone ChIP-seq and ATAC-seq
########################################################################
enrichment_nott_config = dict(
    **{ 'output_dir': 'enrichment/nott',
        'SIGNAL_MANIFEST': '/n/data1/hms/dbmi/park/jluquette/glia/analysis/NOTT.MANIFEST' },
    **enrichment_config
)

module enrichment_nott:
    snakefile: "snakefile.enrichment"
    config: enrichment_nott_config

use rule * from enrichment_nott  as enrichment_nott_*

use rule enrichment_plot from enrichment_nott as enrichment_nott_enrichment_plot with:
    params:
        ignore='enrichment_grid_R_ignore_none',
        group='celltype,mark'
    output:
        expand('plots/enrichment/nott/quantile/{{mutclass}}.{{nquantiles}}quantiles.{output}',
            output=[ 'svg', 'pdf', 'csv' ])
    log:
        'plots/enrichment/nott/quantile/{mutclass}.{nquantiles}quantiles.log'


########################################################################
# Rizzardi et al brain DNA methylation
########################################################################
enrichment_rizzardi_config = dict(
    **{ 'output_dir': 'enrichment/rizzardi',
        'SIGNAL_MANIFEST': '/n/data1/hms/dbmi/park/jluquette/glia/analysis/RIZZARDI.MANIFEST' },
    **enrichment_config
)

module enrichment_rizzardi:
    snakefile: "snakefile.enrichment"
    config: enrichment_rizzardi_config

use rule * from enrichment_rizzardi  as enrichment_rizzardi_*

# Also slightly large; increase memory.
use rule make_qbed_from_bigwig from enrichment_rizzardi as enrichment_rizzardi_make_qbed_from_bigwig with:
    resources:
        mem=15000

use rule enrichment_plot from enrichment_rizzardi as enrichment_rizzardi_enrichment_plot with:
    params:
        ignore='enrichment_grid_R_ignore_none',
        group='region,methyltype,smoothing'
    output:
        expand('plots/enrichment/rizzardi/quantile/{{mutclass}}.{{nquantiles}}quantiles.{output}',
            output=[ 'svg', 'pdf', 'csv' ])
    log:
        'plots/enrichment/rizzardi/quantile/{mutclass}.{nquantiles}quantiles.log'
