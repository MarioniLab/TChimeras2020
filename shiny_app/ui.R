#/bin/R
library(shiny)
library(ggplot2)

load("accessory.RData")

big_plot_width = "1200px"
big_plot_height = "500px"
double_big_plot_height = "1000px"

narrower_plot_width = "650px"

half_plot_width = "600px"
narrower_half_plot_width = "350px"
half_plot_height = "360px"

fluidPage(sidebarLayout(
            sidebarPanel(
              id="sidebar",
              width = 2,
                h3("Plot options"),
                selectizeInput("gene", "Gene", choices = NULL, selected = 26600),
                selectInput(
                  "stage",
                  "Chimera timepoint",
                  choices = c(
                    "E7.5" = "E7.5",
                    "E8.5" = "E8.5"
                  ),
                    selected = "E8.5"
                  ),
                downloadButton("downloadOverview", label = "Chimera expr. vis"),
                downloadButton("downloadSomit", label = "Somitogenesis vis"),
                downloadButton("downloadNMP", label = "NMP vis")
            ),
            mainPanel(
              id = "main",
              width = 10,
              titlePanel(
                "Diverse Routes towards Early Somites in the Mouse Embryo."
              ),
              #put plots here
              tabsetPanel(
                id = "tabs",
                tabPanel(
                  id = "landing",
                  "Landing page",
                  br(),
                  HTML(paste(
                    shiny::tags$b("This is the accompanying interactive server for the paper"),
                    em(
                      "Diverse Routes towards Early Somites in the Mouse Embryo."
                    )
                  )),
                  br(),
                  HTML(
                    "Please note: when you reach this page, the server needs to load some data for your session. We have noticed this is particularly slow on Google Chrome, so please be patient!"
                  ),
                  h4("Tabs:"),
                  HTML(
                    paste(
                      shiny::tags$b("Chimera overview:"),
                      "Reduced dimension plots of the Brachyury chimera data are shown."
                    )
                  ),
                  br(),
                  HTML(
                    paste(
                      shiny::tags$b("Somitogenesis trajectories:"),
                      "Gene expression along the somitogenesis trajectory paths is shown. Specifically, the mean log2 count is shown for all cells in a trajectory at each timepoint."
                    )
                  ),
                  br(),
                  HTML(
                    paste(
                      shiny::tags$b("NMP ordering:"),
                      "Gene expression along the atlas E8.5 NMP cell ordering is shown. Specifically, the mean log2 count is shown for cells in a sliding window across the ordering. There are 100 positions that the window moves across, and each window contains cells from 5% of the trajectory on either side of the window center."
                    )
                  ),
                  h4("Options:"),
                  HTML(paste(
                    shiny::tags$b("Cell subset:"),
                    "Select the chimera timepoints you want to plot (Chimera overview only)."
                  )),
                  br(),
                  HTML(
                    paste(
                      shiny::tags$b("Gene:"),
                      "Select the gene (MGI) to use for expression plots."
                    )
                  ),
                  h4("Other notes:"),
                  HTML(
                    paste(
                      "To report any issues please contact John Marioni at",
                      em("marioni {at} ebi.ac.uk.")
                    )
                  ),
                  br(),
                  HTML(
                    "Use the <a href=\"https://doi.org/doi:10.18129/B9.bioc.MouseGastrulationData\">MouseGastrulationData</a> R package to download the processed data."
                  ),
                ),
                tabPanel(
                  "Chimera overview",
                  id = "overview",
                  plotOutput("overview", width = big_plot_width, height = double_big_plot_height)
                ),
                tabPanel(
                  "Somitogenesis trajectories",
                  id = "somit",
                  plotOutput("somit", width = half_plot_width, height = half_plot_height),
                ),
                tabPanel(
                  "NMP ordering",
                  id = "nmp",
                  plotOutput("nmp", width = half_plot_width, height = half_plot_height),
                )
))))
