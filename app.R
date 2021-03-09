#===============================================================================
# Shiny wep application for an interactive use of 'methylearning' package.
#
# Author: David Piñeyro <davidpvva@hotmail.com>
# Linecese: GPL-3
# Version: 1.0.0
# Date: 2018-05-19
#===============================================================================
library(shiny)
library(methylearning)
library(data.table)
library(GEOquery)
library(Biobase)
library(parallel)
library(xtable)
#library(e1071)
library(shinythemes)

#===============================================================================
# Global variables.
#===============================================================================
# Maximum input file size (in MB).
max_file_size <- 2000
options(shiny.maxRequestSize=max_file_size*1024^2)

#===============================================================================
# Functions.
#===============================================================================
# List fs_methods.
#
# From a character vector with methylearning fs_methods names, this function
# creates a list of these methods, but named with the normal method's name.
# It is created to display regular methods names while keep using methylearnig
# names internally.
#
list_fs_methods <- function(fs_methods) {
  normal_names <- sapply(fs_methods, function(x) {
    if (x == "random_fs") {
      "Random"
    } else if (x == "anova_fs") {
      "Anova"
    } else if (x == "limma_fs") {
      "Limma"
    } else if (x == "correlation_based_fs") {
      "CFS"
    } else if (x == "information_gain_fs") {
      "Information Gain"
    } else if (x == "relief_fs") {
      "Relief"
    } else if (x == "mRMR_fs") {
      "mRMR"
    } else if (x == "ga_fs") {
      "GA"
    } else if (x == "mseq_fs") {
      "mseq"
    } else if (x == "boruta_fs") {
      "boruta"
    }})
  new_list <- as.list(fs_methods)
  names(new_list) <- normal_names
  return(new_list)
}

# List cls_methods.
#
# From a character vector with methylearning cls_methods names, this function
# creates a list of these methods, but named with the normal method's name.
# It is created to display regular methods names while keep using methylearnig
# names internally.
#
list_cls_methods <- function(cls_methods) {
  normal_names <- sapply(cls_methods, function(x) {
    if (x == "knn") {
      "KNN"
    } else if (x == "C5.0") {
      "C5.0"
    } else if (x == "rf") {
      "RandomForest"
    } else if (x == "svmLinear") {
      "svmLinear"
    } else if (x == "svmRadial") {
      "svmRadial"
    } else if (x == "nnet") {
      "nnet"
    } else if (x == "glmnet") {
      "glmnet"
    } else if (x == "lda") {
      "lda"
    }})
  new_list <- as.list(cls_methods)
  names(new_list) <- normal_names
  return(new_list)
}
#===============================================================================
# UI.
#===============================================================================
ui <- fluidPage(
  # Main CSS theme.
  theme = shinytheme("sandstone"),
  # App title
  fluidRow(
    column(11,
           titlePanel(paste0("methylearning: machine learning for ",
                             "DNA methylation (v.1.0.0)"))),
    column(1,
           actionButton(inputId = "exit_button", label = "Exit"))
  ),
  # Add panels. The app is divided into three different panels:
  #   Panel 1: data loading.
  #   Panel 2: feature selection.
  #   Panel 3: classification.
  tabsetPanel(
    #===========================================================================
    # Data Loading
    #
    # In this panel the user specify the input data. There are three
    # possible choices:
    #   User file (csv).
    #   GEO accesion number.
    #   Demo data.
    tabPanel("Data Loading",
      sidebarLayout(
        ### sidebarPanel to select options. =====
        sidebarPanel(
          radioButtons(inputId = "data_input",
                       label = h3("Data uploading options"),
                       choices = list(
                             "None selected" = 1,
                             "CSV file from your system" = 2,
                             "GEO (Gene Expression Omnibus) accession code" = 3,
                             "Demo data set" = 4),
                       selected = 1),
          ### Input: Select a file from your system. =====
          conditionalPanel(
            condition = "input.data_input == 2",
            fileInput(inputId = "user_file",
                       label = h3("Choose csv file from your system ",
                                  "(limit 2 GB)"),
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
            helpText("The input file should be either comma or tab separated,",
                     "with samples as rows and features (array probes or ",
                     "methylation positions) as columns. Values should be ",
                     "beta values for microarray data or methylation ",
                     "ratios for WGBS (whole genome bisulfite sequencing) ",
                     "data. The last column will be taken as the labels for ",
                     "the classification."),
            br(),
            fluidRow(
              column(12,
                     h3("Options:"))
            ),
            checkboxInput(inputId = "sample_names_first_row",
                          label = "Use first column as sample names",
                          value = TRUE),
            checkboxInput(inputId = "user_split_data",
                          label = "Split the data into training and test",
                          value = FALSE),
            conditionalPanel(
              condition = "input.user_split_data",
              textInput(inputId = "user_splitting",
                        label = "Select the proportion of the data to use for training",
                        placeholder = "E.g. 0.66, 2/3, 1/2, 0.5, ..."),
              actionButton(inputId = "confirm_user_split", "Confirm")
            ),
            conditionalPanel(
              condition = "!(input.user_split_data)",
              actionButton(inputId = "confirm_user_no_split", "Confirm")
            )
          ),
          ### Input: Select a GEO accession code. =====
          conditionalPanel(
            condition = "input.data_input == 3",
            textInput(inputId = "geo",
                      label = "Enter a valid GEO accession number",
                      placeholder = "E.g. GSE69229, ..."),
            helpText("VERY IMPORTANT: this method is currently implemented ",
                     "only for data from Illumina Infinium HumanMethylation ",
                     "array (27k, 450k or EPIC). The accesion code should ",
                     "correspond to a single series (one platform)."),
            textInput(inputId = "geo_labels_column",
                      label = "Enter the name of class labels column.",
                      placeholder = paste("Enter a valid column name from ",
                                          "phenoData. E.g. outcome:ch1, ...")),
            helpText("Enter the name of the column to be used as class ",
                     "labels, from", code("phenoData"), "table. This is specific",
                     "for each experiment, depending on the authors choice,",
                     "and may require manual inspection of the 'Series Matrix",
                     "File' available at GEO."),
            checkboxInput(inputId = "geo_split_data",
                          label = "Split the data into training and test",
                          value = FALSE),
            conditionalPanel(
              condition = "input.geo_split_data",
              textInput(inputId = "geo_splitting",
                        label = "Select the proportion of the data to use for training",
                        placeholder = "E.g. 0.66, 2/3, 1/2, 0.5, ..."),
              actionButton(inputId = "confirm_geo_split", "Confirm")
            ),
            conditionalPanel(
              condition = "!(input.geo_split_data)",
              actionButton(inputId = "confirm_geo_no_split", "Confirm")
            )
          ),
          ### Input: Select a demo data set. =====
          conditionalPanel(
            condition = "input.data_input == 4",
            selectInput(inputId = "demo",
                        label = "Select a demo data set",
                        choices = list("Methylation array demo" = 1),
                                      # "WGBS demo" = 2),
                        selected = 1),
            helpText(
              p("This demo data set is a subset of 1000 randomly selected",
                "probes from the publicly available GEO dataset: ",
              a("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69229")),
              p("It was generted using Illumina HumanMethylation array ",
                "450k platform and consists of 52 methylation profiles from ",
                "childhood Acute Lymphoblastic Leukemia (ALL) patients. ",
                "Samples were labeled as 'long.term.remission' and ",
                "'went.on.to.relapse'")),
            checkboxInput(inputId = "demo_split_data",
                          label = "Split the data into training and test",
                          value = FALSE),
            conditionalPanel(
              condition = "input.demo_split_data",
              textInput(inputId = "demo_splitting",
                        label = "Select the proportion of the data to use for training",
                        placeholder = "E.g. 0.66, 2/3, 1/2, 0.5, ..."),
              actionButton(inputId = "confirm_demo_split", "Confirm")
            ),
            conditionalPanel(
              condition = "!(input.demo_split_data)",
              actionButton(inputId = "confirm_demo_no_split", "Confirm")
            )
          ),
          hr(),
          h3("Global Random Seed"),
          helpText("Select a random seed to generate reproducible results."),
          # Select one global random seed to reproduce results
          numericInput(inputId = "global_random_seed",
                       label = "", value = 1234, min = 1, max = 99999999,
                       width = '50%')
        ),
        ### Main panel for displaying outputs. =====
        mainPanel(
          fluidRow(
            column(12,
                   h3("Summary of input data"),
                   verbatimTextOutput("input_summary")
                  )),
          fluidRow(
            column(12,
                   h3("First rows and columns of the input data"),
                   tableOutput("sample_table")
                  ))
        )
      )
    ),
    #===========================================================================
    # Feature Selection
    #
    # In this panel input_data generated in the previous panel will be used
    # for feature selection, depending on user selected options.
    tabPanel("Feature Selection",
      sidebarLayout(
        ### sidebarPanel to select options. =====
        sidebarPanel(
          fluidRow(
            column(12,
                   h2("Feature Selection"),
                   helpText("Please, select one or more Feature Selection",
                            "methods to be applied to your data set.")
            )
          ),
          fluidRow(
            column(6,
                   checkboxGroupInput(inputId = "fs_methods",
                                      label = "Available FS Methods",
                                      choices = list("Random" = 1,
                                                     "Anova" = 2,
                                                     "Limma" = 3,
                                                     "CFS" = 4,
                                                     "Information Gain" = 5,
                                                     "Relief" = 6,
                                                     "mRMR" = 7,
                                                     "GA" = 8,
                                                     "mseq" = 9,
                                                     "boruta" = 10),
                                      selected = character(0)
                   )
            ),
            column(6,
                   numericInput(inputId = "selection_size",
                                label = "Selection Size",
                                value = 30,
                                min = 2),
                   helpText("Maximum number of features to be selected."),
                   hr(),
                   numericInput(inputId = "cores",
                                label = "Number of Processors to Use",
                                value = 1,
                                min = 1,
                                max = detectCores()),
                   helpText("Number of CPU processors to be used in parallel",
                            "computations.")
            )
          ),
          actionButton(inputId = "confirm_fs_method", "Perform FS"),
          helpText("It may take a while, depending on which methods ",
                   "were selected."),
          hr(),
          uiOutput(outputId = "fs_results_control")

        ),
        ### Main panel for displaying outputs. =====
        mainPanel(
          verbatimTextOutput("fs_results_text"),
          dataTableOutput("fs_results_table"),
          verbatimTextOutput("fs_result_plot_error"),
          plotOutput("fs_result_plot")
        )
      )
    ),
    #===========================================================================
    # Classification
    #
    # In this panel fs generated in the previous panel will be used for
    # classification, depending on user selected options.
    tabPanel("Classification",
      sidebarLayout(
        ### sidebarPanel to select options. =====
        sidebarPanel(
          fluidRow(
            column(12,
                   h2("Classification")
            )
          ),
          # The initial UI is rendered depending on the fs_methods used. In this
          # way, only applicable methods appear. See server.
          uiOutput(outputId = "cls_selection"),
          # Next, the non conditional part is rendered.
          hr(),
          h3("Options"),
          fluidRow(
            column(4,
                   numericInput(inputId = "cv_folds",
                                label = "Cross-validation folds",
                                value = 10,
                                min = 2,
                                max = 100)
            ),
            column(4,
                   numericInput(inputId = "repeats",
                                label = "Repetitions",
                                value = 3,
                                min = 1,
                                max = 100)
            ),
            column(4,
                   numericInput(inputId = "tune_length",
                                label = "Tune length",
                                value = 10,
                                min = 1,
                                max = 50)
            )
          ),
          fluidRow(
            column(4,
                   selectInput(inputId = "metric",
                               label = "Select a ranking metric",
                               choices = list("Kappa" = "Kappa",
                                              "Accuracy" = "Accuracy"),
                               selected = 1)
            ),
            column(4,
                   checkboxInput(inputId = "test_eval",
                                 label = paste("Evaluate performance on test ",
                                               "data (if data was previously ",
                                               "partitioned)"),
                                 value = FALSE)
            ),
            column(4,
                   numericInput(inputId = "cores_cls",
                                label = "Number of Processors to Use",
                                value = 1,
                                min = 1,
                                max = detectCores())
            )
          ),
          actionButton(inputId = "confirm_cls_method","Perform Classification"),
          helpText("It may take a while, depending on which methods ",
                   "were selected."),
          # New UI to control results exploration. See server.
          uiOutput(outputId = "cls_results_control")
        ),
        ### Main panel for displaying outputs. =====
        mainPanel(
          verbatimTextOutput("cls_results_text"),
          dataTableOutput("cls_results_table_acc"),
          dataTableOutput("cls_results_table_kpp"),
          plotOutput("cls_results_plot")
        )
      )
    )
  )
)

#===============================================================================
# Server Logic.
#===============================================================================
server <- function(input, output) {
  #=============================================================================
  # Display the welcome message.
  showModal(modalDialog(
    size = "l",
    footer = modalButton("Start"),
    h1("Welcome to methylearning Shiny app!", align = "center"),
    h2("methylearning R package", style = "color:MidnightBlue"),
    p("This is the companion", em("Shiny"), "app for", code("methylearning"),
      "package. It is designed to be a complete framework for",
      "machine learning studies using DNA methylation data."),
    p("The current implementation facilitates the application of several",
      "feature selection (FS) algorithms and several supervised",
      "classification algorithms, providing a wealth of numeric and graphic",
      "performance metrics to be used in the process of model creation."),
    h2("Usage", style = "color:MidnightBlue"),
    p("To use this app you have to navigate the panels following the order,",
      "from left to right:", strong("DATA LOADING"), "|",
      strong("FEATURE SELECTION"), "|", strong("CLASSIFICATION")),
    h3("Data Loading", style = "color:DeepSkyBlue"),
    p("You can choose either to upload your own data, download it from ",
      "Gene Expression Omnibus (GEO) database (only for DNA methylation ",
      "arrays), or use the demo datasaet."),
    p("Once your data is uploaded and you see the summary information, you",
      "can proceed with the next panel."),
    h3("Feature Selection", style = "color:DeepSkyBlue"),
    p("In this panel you can choose the feature selection algorithms to ",
      "apply to your data set, the maximum numbers of features to be selected and",
      "the CPU cores to be used for parallel computations."),
    p("Once feature selection is completed, a bunch of new options will be",
      "displayed, allowing you to study the results of your feature selection and",
      "also store them using the 'Download' button."),
    h3("Classification", style = "color:DeepSkyBlue"),
    p("Finally, once feature selection is completed, the available options for",
      "classification will be displayed."),
    p("After classification completion, a new set of options will be displayed,",
      "allowing you to check the performance for each of the feature selection",
      "+ classification combinations. Several graphical and numerical measures",
      "are available to display and download. The best performant model in",
      "terms of 'accuracy' or 'Kappa' statistic can be easily retrieved to be",
      "used in an", code("R"), "programming language enviroment."),
    p(strong("NOTE:"), "please, be patient as some of the computations can",
      "take a long to complete, depending on your data set and the options",
      "selected.")
  ))
  #=============================================================================
  # Data loading
  #
  # This section contains the panel 1 server logic.
  # First, observe whether exit button was pressed.
  observe({
    if(input$exit_button > 0){
      stopApp("App stopped!")
    }
  })
  # Set the global random seed.
  my_seed <- NULL
  makeReactiveBinding("my_seed")
  observeEvent(input$global_random_seed, {
    my_seed <<- input$global_random_seed
  })

  # Generate a reactive placeholder to store 'ml_data' object, depending
  # on input paramenters.
  input_data <- NULL
  makeReactiveBinding("input_data")
  ### User input data processing for panel 1. =====
  observeEvent(input$confirm_user_split, {
    req(input$user_file)
    req(input$user_splitting)
    data_splitting <- eval(parse(text = input$user_splitting))
    tryCatch(
      {
        df <- fread(input = input$user_file$datapath, stringsAsFactors = TRUE,
                    verbose =FALSE, showProgress = FALSE, data.table = FALSE)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs.
        stop(safeError(e))
      }
    )
    # Sample names.
    if (input$sample_names_first_row) {
      rownames(df) <- df[, 1]
      df <- df[, -1]
    }
    # Generate ml_data object. Assuming last column as class labels
    ml <- ml_data(df, labels_column = dim(df)[2])
    # Data splitting.
    set.seed(isolate({my_seed}))
    ml <- split_data(ml, data_splitting)
    input_data <<- ml
  })

  observeEvent(input$confirm_user_no_split, {
    req(input$user_file)
    tryCatch(
      {
        df <- fread(input = input$user_file$datapath, stringsAsFactors = TRUE,
                    verbose =FALSE, showProgress = FALSE, data.table = FALSE)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs.
        stop(safeError(e))
      }
    )
    # Sample names.
    if (input$sample_names_first_row) {
      rownames(df) <- df[, 1]
      df <- df[, -1]
    }
    # Generate ml_data object. Assuming last column as class labels
    ml <- ml_data(df, labels_column = dim(df)[2])
    input_data <<- ml
  })

  observeEvent(input$confirm_geo_split, {
    req(input$geo)
    req(input$geo_labels_column)
    req(input$geo_splitting)
    data_splitting <- eval(parse(text = input$geo_splitting))
    # Uploading GEO data.
    withProgress(message = 'Downloading from GEO', value = 0, {
      incProgress(0.3)
      gse <- getGEO(input$geo)
      incProgress(0.3)
      ml <- get_GEO_methylarray(gse = gse,
                                target_name = input$geo_labels_column)
      incProgress(0.4)
    })
    # Data splitting.
    set.seed(isolate({my_seed}))
    ml <- split_data(ml, data_splitting)
    input_data <<- ml
  })

  observeEvent(input$confirm_geo_no_split, {
    req(input$geo)
    req(input$geo_labels_column)
    # Uploading GEO data. 
    withProgress(message = 'Downloading from GEO', value = 0, {
      incProgress(0.3)
      gse <- getGEO(input$geo)
      incProgress(0.3)
      ml <- get_GEO_methylarray(gse = gse,
                                target_name = input$geo_labels_column)
      incProgress(0.4)
    })
    input_data <<- ml
  })

  observeEvent(input$confirm_demo_split, {
    req(input$demo)
    req(input$demo_splitting)
    data_splitting <- eval(parse(text = input$demo_splitting))
    if (input$demo == 1) {
      demo_file <- "demo_data/demo_1.RData"
    } else if (input$demo == 2) {
      demo_file <- "demo_data/demo_2.RData"
    }
    tryCatch(
      {
        load(demo_file)
        ml <- ml_data(df, labels_column = dim(df)[2])
      },
      error = function(e) {
        # return a safeError if a parsing error occurs.
        stop(safeError(e))
      }
    )
    # Data splitting.
    set.seed(isolate({my_seed}))
    ml <- split_data(ml, data_splitting)
    input_data <<- ml
  })

  observeEvent(input$confirm_demo_no_split, {
    req(input$demo)
    if (input$demo == 1) {
      demo_file <- "demo_data/demo_1.RData"
    } else if (input$demo == 2) {
      demo_file <- "demo_data/demo_2.RData"
    }
    tryCatch(
      {
        load(demo_file)
        ml <- ml_data(df, labels_column = dim(df)[2])
      },
      error = function(e) {
        # return a safeError if a parsing error occurs.
        stop(safeError(e))
      }
    )
    input_data <<- ml
  })

  ### Panel 1 reactivity. =====
  observe({
    output$input_summary <- renderPrint({
      if (is.null(input_data)) {
        cat("Input data has not been loaded yet!")
      } else {
        summary(input_data)
      }
    })
  })
  observe({
    output$sample_table <- renderTable({
    input_data$training_data[1:4, 1:8]
    })
  })

  #=============================================================================
  # Feature Selection Part 1: compute FS.
  #
  # This section contains the panel 2, part 1 server logic.
  # Initialize this variable as NULL to know whether it was really initialized.
  fs <- NULL
  fs <- eventReactive(input$confirm_fs_method, {
    if (is.null(input_data)) {
      return(cat("Input data has not been loaded yet!\n"))
    }
    # Assert data was entered.
    req(input$fs_methods)
    req(input$selection_size)
    req(input$cores)
    fs_methods <- available_fs_methods()$all[as.numeric(input$fs_methods)]
    # Perform Feature Selection.
    set.seed(isolate({my_seed}))
    # Create ml_fs object with added progress bar.
    withProgress(message = 'Performing Feature Selection', value = 0, {
      incProgress(0.3)
      f <- ml_fs(ml_data_obj = input_data, 
                 fs_methods = fs_methods, 
                 selection_size = input$selection_size, 
                 cores = input$cores)
      incProgress(0.7)
    })

    return(f)
  })
  # Render a new UI panel after FS computation to select graph options.
  output$fs_results_control <- renderUI({
    # Render a new UI based on selection methods used. This new UI will be in
    # charge of render graphics.
    fs_methods <- list_fs_methods(selection_methods(fs()))
    fluidRow(
      column(12,
             h3("Feature Selection Results Exploration"),
             radioButtons(inputId = "fs_choice",
                          label = "Select one of the following actions:",
                          choices = list("Show FS summary" = 1,
                                         "Show FS results" = 2,
                                         "Show FS summary table" = 3,
                                         "Plot Venn diagram" = 4,
                                         "Plot computation time" = 5),
                          selected = 1),
             conditionalPanel(
               condition = "input.fs_choice == 4",
               checkboxGroupInput(inputId = "fs2venn",
                                  label = paste("Select FS methods results to ",
                                                "compare (from 2 up to 5)"),
                                  choices = fs_methods)
             ),
             fluidRow(
               column(6,
                      actionButton(inputId = "confirm_fs_choice",
                                   label = "Confirm")
               ),
               column(6,
                      downloadButton("download_fs_results", "Download")
               )
             )
      )
    )
  })
  #=============================================================================
  # Feature Selection Part 2: explore/download FS results.
  #
  # This section contains the panel 2, part 2 server logic.
  fs_exploration <- eventReactive(input$confirm_fs_choice, {
    # This event reactive reacts to confirm_fs_choice actionButton and
    # returns a different output, depending on input$fs_choice.
    req(input$fs_choice)
    if (input$fs_choice == 1) {
      return(1)
    } else if (input$fs_choice == 2) {
      return(2)
    } else if (input$fs_choice == 3) {
      return(3)
    } else if (input$fs_choice == 4) {
      req(input$fs2venn)
      return(4)
    } else if (input$fs_choice == 5) {
      return(5)
    }
  })

  output$fs_results_text <- renderPrint({
    fs_choice <- fs_exploration()
    if (fs_choice == 1) {
      summary(fs())
    } else if(fs_choice == 2) {
      selection_results(fs())
    }
  })

  output$fs_results_table <- renderDataTable({
    fs_choice <- fs_exploration()
    if (fs_choice == 3) {
      xtable(cbind(Sample = rownames(selection_df(fs())), selection_df(fs())))
    }
  })

  output$fs_result_plot_error <- renderPrint({
    fs_choice <- fs_exploration()
    if(fs_choice == 4) {
      fs2plot <- isolate({input$fs2venn})
      l_fs <- length(fs2plot)
      if (l_fs < 2 | l_fs > 5) {
        return(cat("Select from 2 up to 5 FS methods to plot!"))
      }
    }
  })

  output$fs_result_plot <- renderPlot({
    fs_choice <- fs_exploration()
    if(fs_choice == 4) {
      fs2plot <- isolate({input$fs2venn})
      l_fs <- length(fs2plot)
      if (l_fs > 1 & l_fs < 6) {
        plot_venn_fs(fs(), fs_methods = fs2plot,
                     fill = rainbow(length(fs2plot)),
                     category = fs2plot)
      }
    } else if (fs_choice == 5) {
    plot_fs_time(fs())
    }
  })

  # Download content. Summary cannot be downloaded (but could be directly
  # copied).
  output$download_fs_results <- downloadHandler(
    filename = function(){
      if (fs_exploration() == 4 | fs_exploration() == 5) {
        paste0("fs_plot_", Sys.Date(), ".png")
      } else {
        paste0("fs_results_", Sys.Date(), ".csv")
      }
    },
    content = function(file) {
      if (fs_exploration() == 1 | fs_exploration() == 2) {
        # Correct for different lengths of FS.
        sel_fs <- selection_results(fs())
        lengths_vec <- unlist(lapply(sel_fs, length))
        if (length(unique(lengths_vec)) > 1) {
          # Add NAs.
          max_features <- max(lengths_vec)
          to_add_NAs <- sapply(lengths_vec, function(x) max_features - x)
          for (i in 1:length(sel_fs)) {
            sel_fs[[i]] <- c(sel_fs[[i]], rep(NA, as.numeric(to_add_NAs[i])))
          }
        }
        write.csv(sel_fs, file)
      } else if (fs_exploration() == 3) {
        fwrite(selection_df(fs()), file, row.names = TRUE)
      } else if (fs_exploration() == 4) {
        fs2plot <- isolate({input$fs2venn})
        l_fs <- length(fs2plot)
        if (l_fs > 1 & l_fs < 6) {
          png(file)
          plot_venn_fs(fs(), fs_methods = fs2plot,
                       fill = rainbow(length(fs2plot)),
                       category = fs2plot)
          dev.off()
        }
      } else if (fs_exploration() == 5) {
        png(file)
        pl <- plot_fs_time(fs())
        print(pl)
        dev.off()
      }
    }
  )

  #=============================================================================
  # Classification Part 1: perform classificaiton
  #
  # This section contains the panel 3 server logic.
  # First, the methods selection part of the panel should be displayed when
  # fs_methods become available.
  output$cls_selection <- renderUI({
    fs_methods <- list_fs_methods(selection_methods(fs()))
    fluidRow(
      column(12,
             helpText(
               "Please, select a combination of obtained feature selections and",
               "classification methods to be applied to each feature set.")
      ),
      column(6,
             checkboxGroupInput(inputId = "fs_methods_performed",
                                label = "Feature Sets Available",
                                choices = fs_methods,
                                selected = character(0)
             )
      ),
      column(6,
             checkboxGroupInput(inputId = "cls_methods",
                                label = "Available Classification Methods",
                                choices = list("KNN" = "knn",
                                               "C5.0" = "C5.0",
                                               "RandomForest" = "rf",
                                               "svmLinear" = "svmLinear",
                                               "svmRadial" = "svmRadial",
                                               "nnet" = "nnet",
                                               "glmnet" = "glmnet",
                                               "lda" = "lda"),
                                selected = character(0)
             )
      )
    )
  })

  cls <- eventReactive(input$confirm_cls_method, {
    if (is.null(input_data) | is.null(fs)) {
      return(cat("Input data has not been loaded yet!\n"))
    }
    # Assert data was entered.
    req(input$fs_methods_performed)
    req(input$cls_methods)
    req(input$cv_folds)
    req(input$repeats)
    req(input$tune_length)
    req(input$metric)
    req(input$cores_cls)
    # Perform classification
    set.seed(isolate({my_seed}))
    # Use a progress bar indicator.
    withProgress(message = 'Performing Classification', value = 0, {
      incProgress(0.3)
      c <- ml_cls(ml_fs_obj = fs(), 
                  cls_methods = input$cls_methods, 
                  fs_methods = input$fs_methods_performed,
                  cv_folds = input$cv_folds, 
                  rep = input$repeats, 
                  tune_length = input$tune_length,
                  metric = input$metric, 
                  test_eval = input$test_eval, 
                  cores = input$cores_cls)
      incProgress(0.7)
    })
    return(c)
  })

  # Render a new UI panel after CLS computation to select graph options.
  output$cls_results_control <- renderUI({
    # Render a new UI based on selection methods used. This new UI will be in
    # charge of render graphics.
    fs_methods <- list_fs_methods(cls()$fs_methods)
    cls_methods <- list_cls_methods(cls()$cls_methods)
    evaluation_performed <- ifelse(class(cls()$test_evaluation) == "list",
                                   TRUE,
                                   FALSE)
    if (evaluation_performed) {
      fluidRow(
        column(12,
               h3("Feature Selection + Classification Results"),
               radioButtons(inputId = "cls_choice",
                            label = "Select one of the following actions:",
                            choices = list(
                          "Show classification summary" = 1,
                          "Show training results" = 2,
                          "Show test results" = 3,
                          "Show best model from training data" = 4,
                          "Show best feature set from training data" = 5,
                          "Show best model from test data" = 6,
                          "Show best feature set from test data" = 7,
                          "Show confusion matrix for test evaluation" = 8,
                          "Plot training validation accuracy" = 9,
                          "Plot training validation Kappa" = 10,
                          "Plot test accuracy" = 11,
                          "Plot test Kappa" = 12,
                          "Plot ROC curves (only for 2 classes datasets)" = 13,
                          "Plot computation time" = 14),
                              selected = 1
               ),
               conditionalPanel(
                 condition = "input.cls_choice == 8",
                 fluidRow(
                   column(6,
                          radioButtons(inputId = "fs2confusionMatrix",
                                       label = paste("Select one FS method"),
                                       choices = fs_methods)
                   ),
                   column(6,
                          radioButtons(inputId = "cls2confusionMatrix",
                                       label = paste("Select one classification method"),
                                       choices = cls_methods)
                   )
                 )
               ),
               conditionalPanel(
                 condition = "input.cls_choice == 13",
                 fluidRow(
                   column(6,
                          radioButtons(inputId = "fs2roc",
                                       label = paste("Select one FS method"),
                                       choices = fs_methods)
                   ),
                   column(6,
                          radioButtons(inputId = "cls2roc",
                                       label = paste("Select one classification method"),
                                       choices = cls_methods)
                   )
                 )
               ),
               conditionalPanel(
                 condition = "input.cls_choice == 14",
                 fluidRow(
                   column(12,
                          radioButtons(inputId = "time2plot",
                                       label = paste("Select one option"),
                                       choices = list(
                      "Computation time for FS + classification" = "sum",
                      "Computation time only for classification" = "cls_only")
                                       )
                   )
                 )
               ),
               fluidRow(
                 column(6,
                        actionButton(inputId = "confirm_cls_choice",
                                     label = "Confirm")
                 ),
                 column(6,
                        downloadButton("download_cls_results",
                                       "Download Displayed Results")
                 )
               ),
               hr(),
               fluidRow(
                 column(12,
                        downloadButton(
                          "download_best_model",
                          "Download Best Model From Training Validation")
                 )
               )
        )
      )

    } else {
      fluidRow(
        column(12,
               h3("Feature Selection + Classification Results"),
               radioButtons(inputId = "cls_choice",
                            label = "Select one of the following actions:",
                            choices = list(
                          "Show classification summary" = 1,
                          "Show training results" = 2,
                          "Show best model from training data" = 4,
                          "Show best feature set from training data" = 5,
                          "Plot training validation accuracy" = 9,
                          "Plot training validation Kappa" = 10,
                          "Plot ROC curves (only for 2 classes datasets)" = 13,
                          "Plot computation time" = 14),
                            selected = 1
               ),
               conditionalPanel(
                 condition = "input.cls_choice == 13",
                 fluidRow(
                   column(6,
                          radioButtons(inputId = "fs2roc",
                                       label = paste("Select one FS method"),
                                       choices = fs_methods)
                   ),
                   column(6,
                          radioButtons(inputId = "cls2roc",
                                       label = paste("Select one classification method"),
                                       choices = cls_methods)
                   )
                 )
               ),
               conditionalPanel(
                 condition = "input.cls_choice == 14",
                 fluidRow(
                   column(12,
                          radioButtons(inputId = "time2plot",
                                       label = paste("Select one option"),
                                       choices = list(
                      "Computation time for FS + classification" = "sum",
                      "Computation time only for classification" = "cls_only")
                                       )
                   )
                 )
               ),
               fluidRow(
                 column(6,
                        actionButton(inputId = "confirm_cls_choice",
                                     label = "Confirm")
                 ),
                 column(6,
                        downloadButton("download_cls_results",
                                       "Download Displayed Results")
                 )
               ),
               hr(),
               fluidRow(
                 column(12,
                        downloadButton(
                          "download_best_model",
                          "Download Best Model From Training Validation")
                 )
               )
        )
      )
    }
  })
  cls_exploration <- eventReactive(input$confirm_cls_choice, {
    # This event reactive reacts to confirm_cls_choice actionButton and
    # returns a different output, depending on input$cls_choice.
    req(input$cls_choice)
    # Ensure that conditional values are also there, when necessary.
    if (input$cls_choice == 8) {
      req(input$fs2confusionMatrix)
      req(input$cls2confusionMatrix)
      return(8)
    } else if (input$cls_choice == 13) {
      req(input$fs2roc)
      req(input$cls2roc)
      return(13)
    } else if (input$cls_choice == 14) {
      req(input$time2plot)
      return(14)
    } else {
      return(input$cls_choice)
    }
  })

  # Generate the outputs.
  output$cls_results_text <- renderPrint({
    # Outputs all the text renders.
    if (cls_exploration() == 1) {
      summary(cls())
    } else if (cls_exploration() == 4) {
      get_best_model(cls())
    } else if (cls_exploration() == 5) {
      get_best_feature_set(cls())
    } else if (cls_exploration() == 6) {
      get_best_model_test(cls())
    } else if (cls_exploration() == 7) {
      get_best_feature_set_test(cls())
    } else if (cls_exploration() == 8) {
      get_confusionMatrix_test(
        cls(),
        fs_method = isolate({input$fs2confusionMatrix}),
        cls_method = isolate({input$cls2confusionMatrix}))
    }
  })
  output$cls_results_table_acc <- renderDataTable({
    # Renders Data Tables.
    if (cls_exploration() == 2) {
      xtable(cbind(FSmethod.Accuracy = rownames(get_training_results(cls())[[1]]),
                   get_training_results(cls())[[1]]))
    } else if (cls_exploration() == 3) {
      xtable(cbind(FSmethod.Accuracy = rownames(get_test_results(cls())[[1]]), 
                   get_test_results(cls())[[1]]))
    }
  })
  output$cls_results_table_kpp <- renderDataTable({
    # Renders Data Tables.
    if (cls_exploration() == 2) {
      xtable(cbind(FSmethod.Kappa = rownames(get_training_results(cls())[[2]]), 
                   get_training_results(cls())[[2]]))
    } else if (cls_exploration() == 3) {
      xtable(cbind(FSmethod.Kappa = rownames(get_test_results(cls())[[2]]), 
                   get_test_results(cls())[[2]]))
    }
  })
  output$cls_results_plot <- renderPlot({
    # Renders all the plots.
    if (cls_exploration() == 9) {
      plot(cls())
    } else if (cls_exploration() == 10) {
      plot(cls(), mode = "Kappa")
    } else if (cls_exploration() == 11) {
      plot_test(cls(), mode = "Accuracy")
    } else if (cls_exploration() == 12) {
      plot_test(cls(), mode = "Kappa")
    } else if (cls_exploration() == 13) {
      plot(cls(),
           mode = "ROC",
           fs_ROC = isolate({input$fs2roc}),
           cls_ROC = isolate({input$cls2roc}))
    } else if (cls_exploration() == 14) {
      plot(cls(),
           mode = "Time",
           time_mode = isolate({input$time2plot}))
    }
  })
  # Download button behaviour.
  output$download_cls_results <- downloadHandler(
    filename = function(){
      if (cls_exploration() %in% c(9, 10, 11, 12, 13, 14)) {
        paste0("cls_plot_", Sys.Date(), ".png")
      } else {
        paste0("cls_results_", Sys.Date(), ".csv")
      }
    },
    content = function(file) {
      # To produce download content, only selections 2, 3, 5, 7, 9, 10, 11, 12
      # and 13 are taken into account.
      if (cls_exploration() %in% c(1, 2, 4, 6, 8)) {
        # All of them will react as if cls_exploration() == 2.
        # CSVs
        write.csv(get_training_results(cls()), file)
      } else if (cls_exploration() == 3 &
                 class(cls()$test_evaluation) == "list") {
        write.csv(get_test_results(cls()), file)
      } else if (cls_exploration() == 5) {
        write.csv(get_best_feature_set(cls()), file)
      } else if (cls_exploration() == 7 &
                 class(cls()$test_evaluation) == "list") {
        write.csv(get_best_feature_set_test(cls()), file)
      } else if (cls_exploration() == 9) {
        # PLOTS.
        png(file)
        pl <- plot(cls())
        print(pl)
        dev.off()
      } else if (cls_exploration() == 10) {
        # PLOTS.
        png(file)
        pl <- plot(cls(), mode = "Kappa")
        print(pl)
        dev.off()
      } else if (cls_exploration() == 11) {
        # PLOTS.
        png(file)
        pl <- plot_test(cls(), mode = "Accuracy")
        print(pl)
        dev.off()
      } else if (cls_exploration() == 12) {
        # PLOTS.
        png(file)
        pl <- plot_test(cls(), mode = "Kappa")
        print(pl)
        dev.off()
      } else if (cls_exploration() == 13) {
        # PLOTS.
        png(file)
        pl <- plot(cls(),
                   mode = "ROC",
                   fs_ROC = isolate({input$fs2roc}),
                   cls_ROC = isolate({input$cls2roc}))
        print(pl)
        dev.off()
      } else if (cls_exploration() == 14) {
        # PLOTS.
        png(file)
        pl <- plot(cls(),
                   mode = "Time",
                   time_mode = isolate({input$time2plot}))
        print(pl)
        dev.off()
      }
    }
  )
  # Download Best Model
  output$download_best_model <- downloadHandler(
    filename <- paste0("best_model_", Sys.Date(),".RData"),
    content = function(file) {
      best_model <- get_best_model(cls())
      save(best_model, file = file)
    }
  )
}

#===============================================================================
# Create Shiny app
#===============================================================================
shinyApp(ui, server)
