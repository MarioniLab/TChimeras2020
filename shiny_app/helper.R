##This script contains useful data: palettes, dimensions, etc.

## DIMENSIONS
big_plot_width = 9 * 2
big_plot_height = 5 * 2
small_plot_width = 6.5 * 1.5
small_plot_height = 5 * 1.5 * 6.5/9

## PALETTES
#celltypes split according to the logic of grouping
celltype_palette = c("Epiblast" = "#635547",
                     "Primitive Streak" = "#DABE99",
                     "Caudal epiblast" = "#9e6762",

                     "PGC" = "#FACB12",
                     
                     "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     
                     "Visceral endoderm" = "#F6BFCB",
                     "ExE endoderm" = "#7F6874",
                     "ExE ectoderm" = "#989898",
                     "Parietal endoderm" = "#1A1A1A"
                     
)

nmp_palette = c(
       "Spinal cord" = "#5281A8",
       "NMP" = "#5E7958",
       "Caudal Mesoderm" = "#75A430",
       "Presomitic mesoderm" = "#ADD445",
       "Posterior-most somites" = "#DDC333",
       "Dermomyotome" = "#CA681E",
       "Head mesoderm" = "#704919"
       )
nmp_labels = c(
       "Spinal cord" = "Spinal cord",
       "NMP" = "NMP",
       "Caudal Mesoderm" = "Caudal mesoderm",
       "Presomitic mesoderm" = "Presomitic mesoderm",
       "Posterior-most somites" = "Posterior somitic tissues",
       "Dermomyotome" = "Dermomyotome",
       "Head mesoderm" = "Cranial mesoderm"
)