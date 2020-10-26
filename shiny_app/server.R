#!/bin/R
library(shiny)
library(ggplot2)
library(HDF5Array)
library(reshape2)
library(cowplot)
theme_set(theme_cowplot())


#load palettes etc.
source("helper.R")

load("accessory.RData")

count_link = HDF5Array(file = "counts.hdf5", name = "logcounts")
nmp_ordering_link = HDF5Array(file = "nmp_ordering.hdf5", name = "windows")

ant_link = HDF5Array(file = "anterior.hdf5", name = "mean")
post_link = HDF5Array(file = "posterior.hdf5", name = "mean")
nmp_link = HDF5Array(file = "nmp.hdf5", name = "mean")

nmp_trace = readRDS("nmp_trace.rds")
nmp_trace$start100 = nmp_trace$start/max(nmp_trace$end) * 100
nmp_trace$end100 = nmp_trace$end/max(nmp_trace$end) * 100
shinyServer(
  function(input, output, session){

    updateSelectizeInput(session = session, inputId = 'gene', choices = genes[,2], server = TRUE, selected = "T") #DONT remove T, appears to be a bug that it vanishes

    ####
    # Expression UMAPs
    ####
    
    get_counts = reactive({
        counts = as.numeric(count_link[,match(input$gene, genes[,2])])
        counts[meta$stage == input$stage]
    })
    get_coords = reactive({
        umaps[[input$stage]]
    })
    get_celltype = reactive({
        meta[meta$stage == input$stage, "celltype.mapped"]
    })
    get_tomato = reactive({
        tom = meta[meta$stage == input$stage, "tomato"]
        sapply(as.character(tom), switch, "TRUE" = "Tomato+", "FALSE" = "Tomato-")
    })

    plot_celltype_umap = function(){
        pdf = data.frame(
            coord_x = get_coords()[,1],
            coord_y = get_coords()[,2],
            celltype = get_celltype(),
            tomato = get_tomato()
        )
        ggplot(pdf, aes(x = coord_x, y = coord_y, col = celltype)) +
            geom_point() +
            labs(x = "umapX", y = "umapY") + 
            scale_colour_manual(values = celltype_palette, name = "") +
            theme(
                axis.text = element_blank(), 
                axis.ticks = element_blank()
                )  +
            facet_wrap(~tomato) +
            guides(colour = guide_legend(override.aes = list(size=7, alpha = 1))) +
            ggtitle(input$stage)
    }

    plot_expr_umap = function(){
        pdf = data.frame(
            coord_x = get_coords()[,1],
            coord_y = get_coords()[,2],
            expr = get_counts(),
            tomato = get_tomato()
        )
        pdf = pdf[order(pdf$expr),]
        ggplot(pdf, aes(x = coord_x, y = coord_y, col = expr)) +
            geom_point() +
            labs(x = "umapX", y = "umapY") + 
            scale_color_gradient2(
                name = "log2\nnormalised\ncounts",
                mid = "cornflowerblue",
                low = "gray75",
                high = "black",
                midpoint = max(pdf$expr)/2) +
            theme(
                axis.text = element_blank(), 
                axis.ticks = element_blank()
                ) +
            facet_wrap(~tomato) +
            ggtitle(input$gene)
    }

    make_overview_grid = function(){
        plot_grid(plot_celltype_umap(), plot_expr_umap(),
                  ncol = 1, align = "hv", axis = "tblr")
    }
    
    output$overview = renderPlot({
        validate(
        need(input$gene %in% genes[,2],
             "Please select a gene; if you have already selected one, this gene is not in our annotation." )
        )
        make_overview_grid()
    })

    output$downloadOverview <- downloadHandler(
      filename = function() { paste0("Tchimera_overview_", input$stage, "-cells_", input$gene, ".pdf") },
      content = function(file) {
        pdf(file = NULL)
        ggsave(file, plot = make_overview_grid(), device = "pdf", width = big_plot_width, height = big_plot_height)
        dev.off()
      }
    )

    ####
    # Somitogenesis plots
    ####
    stages = c("E6.5", "E6.75", "E7.0", "E7.25", "E7.5",
               "E7.75", "E8.0", "E8.25", "E8.5")
    get_somite_data = reactive({
        index = match(input$gene, genes[,2])
        do.call(rbind, list(
            data.frame(expr = ant_link[,index], ct = "Anterior somites", stage = stages),
            data.frame(expr = post_link[,index], ct = "Posterior somites", stage = stages),
            data.frame(expr = nmp_link[,index], ct = "NMP", stage = stages)))
    })
    plot_somite_data = function(){
        pdf = get_somite_data()
        ggplot(pdf, aes(x = stage, y=expr, group = ct, col = ct)) +
            geom_line(lwd = 2) +
            labs(y = sprintf("Mean %s log2 normalised count", input$gene)) +
            theme(axis.title.x = element_blank()) +
            scale_colour_manual(values = c(
                "Anterior somites" = "#E33C25",
                "Posterior somites" = "#DDC333",
                "NMP" = "#5E7958"
            ), name = "") +
            ggtitle(input$gene)
    }
    output$somit = renderPlot({
        validate(
        need(input$gene %in% genes[,2],
             "Please select a gene; if you have already selected one, this gene is not in our annotation." )
        )
        plot_somite_data()})
    output$downloadSomit <- downloadHandler(
      filename = function() { paste0("somitogenesis_", input$gene, ".pdf") },
      content = function(file) {
        pdf(file = NULL)
        ggsave(file, plot = plot_somite_data(), device = "pdf", width = small_plot_width, height = small_plot_height)
        dev.off()
      }
    )

    ####
    # NMP plots
    ####
    get_nmp_data = reactive({
        index = match(input$gene, genes[,2])
        as.numeric(nmp_ordering_link[,index])
    })

    plot_nmp_data = function(){
        ggplot() +
            geom_rect(data = nmp_trace, 
                mapping = aes(xmin = start100, xmax = end100, 
                              ymin = 0, ymax = -max(get_nmp_data()) * 0.05,
                              fill = call)) +
            geom_path(mapping = aes(x = seq_along(get_nmp_data()), y = get_nmp_data())) +
            labs(x = "NMP ordering", y = sprintf("Average %s log2 normalised counts", input$gene)) +
            theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
            scale_fill_manual(values = nmp_palette, labels = nmp_labels, name = "Most common celltype\nin window") +
            ggtitle(input$gene)
    }
    
    output$nmp = renderPlot({
        validate(
        need(input$gene %in% genes[,2],
             "Please select a gene; if you have already selected one, this gene is not in our annotation." )
        )
        plot_nmp_data()})
    output$downloadNMP <- downloadHandler(
      filename = function() { paste0("NMP_", input$gene, ".pdf") },
      content = function(file) {
        pdf(file = NULL)
        ggsave(file, plot = plot_nmp_data(), device = "pdf", width = small_plot_width, height = small_plot_height)
        dev.off()
      }
    )

})
