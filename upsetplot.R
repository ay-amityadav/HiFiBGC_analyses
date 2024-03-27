# Creates Upsetplots 

library(ggplot2)
library(ComplexUpset)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_dir <- args[2]
c_parameter <- args[3]

set_size = function(w, h, factor=1.5) {
    s = 1 * factor
    options(
        repr.plot.width=w * s,
        repr.plot.height=h * s,
        repr.plot.res=100 / factor,
        jupyter.plot_mimetypes='image/png',
        jupyter.plot_scale=1
    )
}

df <- read.csv(input_file, sep='\t', header=TRUE)

# Rename column names
names(df)[names(df) == "Partial.Complete"] <- "Complete/Partial"
names(df)[names(df) == "BiGSCAPE_class"] <- "BiG-SCAPE Class"

if (ncol(df) == 7) {
    columns <- c('hicanu', 'metaflye', 'hifiasm.meta')
} else {
    columns <- c('hicanu', 'metaflye', 'hifiasm.meta', 'unmapped_reads')
}

################
# Output as pdf
################

# Upsetplot 
output_file = paste(output_dir, "/", "upsetplot", "_", c_parameter, ".pdf", sep='')
pdf(output_file) 
set_size(8, 3)
upset(
    df,
    columns,
    base_annotations=list(
        'Intersection size'=intersection_size(
            counts=TRUE,
            mapping=aes(fill=`Complete/Partial`),
            bar_number_threshold=1.00
        )
    ),
    set_sizes=(
        upset_set_size()
        + theme(axis.ticks.x=element_line())
    ),
    width_ratio=0.360,
 ) + patchwork::plot_layout(heights=c(1.2, 0.5))
dev.off()

# Upsetplot with percentages 
output_file = paste(output_dir, "/", "upsetplot", "_", c_parameter, "_with_percentages", ".pdf", sep='')
pdf(output_file) 
set_size(8, 3)
upset(
    df,
    columns,
    base_annotations=list(
        'Intersection size'=intersection_size(text_mapping=aes(label=paste0(round(
            !!get_size_mode('exclusive_intersection')/(nrow(df) - 1) * 100
            ), '%')),
            text=list(size=2.9),
            width=0.9,
            mapping=aes(fill=`Complete/Partial`),
            bar_number_threshold=1.00
        )
    ),
    set_sizes=(
        upset_set_size()
        + theme(axis.ticks.x=element_line())
    ),
    width_ratio=0.360,
 ) + patchwork::plot_layout(heights=c(1.2, 0.5))
dev.off()

# Upsetplot with clinker_lower_matrix_mean
output_file = paste(output_dir, "/", "upsetplot", "_", c_parameter, "_clinker_lower_matrix_mean", ".pdf", sep='')
pdf(output_file) 
set_size(8, 3)
upset(
    df,
    columns,
    base_annotations=list(
        'Intersection size'=intersection_size(
            counts=TRUE,
            mapping=aes(fill=`Complete/Partial`),
            bar_number_threshold=1.00
        )
    ),
    annotations=list(
        'clinker_lower_matrix_mean'=list(
            aes=aes(x=intersection, y=clinker_lower_matrix_mean),
            geom=geom_boxplot(na.rm=TRUE)
        )
    ),
    set_sizes=(
        upset_set_size()
        + theme(axis.ticks.x=element_line())
    ),
    width_ratio=0.360,
 ) + patchwork::plot_layout(heights=c(0.7, 1.2, 0.5))
dev.off()

# Upsetplot with BiG-SCAPE classes
output_file = paste(output_dir, "/", "upsetplot", "_", c_parameter, "_bigscape_class", ".pdf", sep='')
pdf(output_file) 
set_size(8, 3)
upset(
    df,
    columns,
    base_annotations=list(
        'Intersection size'=intersection_size(
            counts=TRUE,
            mapping=aes(fill=`Complete/Partial`),
            bar_number_threshold=1.00
        )
    ),
    annotations = list(
        'BiG-Scape class'=(
            ggplot(mapping=aes(fill=`BiG-SCAPE Class`))
            + geom_bar(stat='count', position='fill')
            + scale_y_continuous(labels=scales::percent_format())
            + scale_fill_manual(values=c(
                 'NRPS'='#0072B2', 'Others'='#999999',
                 'PKS-NRP_Hybrids'='#56B4E9', 'PKSI'='#009E73',
                 'PKSother'='#F0E442', 'RiPPs'='#E69F00',
                 'Saccharides'='#D55E00', 'Terpene'='#CC79A7'
            ))
            + ylab('Count percentage')
        )
    ),
    set_sizes=(
        upset_set_size()
        + theme(axis.ticks.x=element_line())
    ),
    width_ratio=0.360,
 ) + patchwork::plot_layout(heights=c(0.76, 1.2, 0.5))
dev.off()

