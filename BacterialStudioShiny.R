library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(dada2)
library(stringr)
library(ggplot2)
library(DT)

ui <-
  fluidPage(titlePanel("NGS & BGD Studio"),
            hr(),
            tabsetPanel(tabPanel("NGSstudio",
                                 fluidRow(
                                   navlistPanel(
                                     widths = c(2, 10),
                                     tabPanel(
                                       "Prepare Reads",
                                       br(),
                                       navbarPage(
                                         "Reads Preparation >>",
                                         tabPanel("MergePairs",
                                                  helpText("* If pair-reads were not merged, 
                                                           you should firstly upload that. 
                                                           File format: .fq or .fastq. "),
                                                  helpText(">> vsearch --fastq_mergepairs fastqfile --reverse fastqfile --fastqout outputfile"),
                                                  helpText("* Merge paired-end sequence reads into one sequence. 
                                                           The forward reads are specified as the argument to this 
                                                           option and the reverse reads are specified with the --reverse option.
                                                           The merged sequences are output to the file(s) specified with 
                                                           the --fastqout options. "),
                                                  hr(),
                                                  fluidRow(
                                                    column(
                                                      4,
                                                      h5(strong("Select the folder of pair-end reads to parse the directory...")),
                                                      shinyDirButton(
                                                        "rawPairReadsPath",
                                                        "Select Folder",
                                                        "Upload"
                                                      )
                                                    )
                                                  ),
                                                  hr(),
                                                  h4("> Input the identification symbol..."),
                                                  fluidRow(
                                                    column(
                                                      4,
                                                      textInput(
                                                        inputId = "difsymbolForward",
                                                        label = "Forward (e.g., _R1.fq, and _Forward.fq.):"
                                                      )
                                                    ),
                                                    column(
                                                      4,
                                                      textInput(
                                                        inputId = "difsymbolReverse",
                                                        label = "Reverse (e.g., _R2.fq, and _Reverse.fq.):"
                                                      )
                                                    )
                                                  ),
                                                  hr(),
                                                  helpText("Next merge uploaded pair-reads..."),
                                                  fluidRow(column(4, progressBar(id = "prgsBar_mergePairs", value = 0, display_pct = T))),
                                                  actionBttn(
                                                    inputId = "mergePairs",
                                                    label = "Merge Pair-reads files...",
                                                    style = "jelly",
                                                    color = "primary"
                                                  )
                                                  ),
                                         tabPanel("QualityFilter",
                                                  helpText("* Quality filtering is critical in reducing the abundance 
                                                           and impact of spurious sequences. There is an intrinsic 
                                                           error rate to all sequencing technologies (and polymerases) 
                                                           that will consistently be generating some portion of 
                                                           sequences that vary from their true biological origin, 
                                                           and this can substantially inflate metrics such as 
                                                           richness and diversity. Quality filtering is one of 
                                                           the steps in place to minimize that problem. 
                                                           File format: .fq or .fastq.  or .fasta"),
                                                  helpText(">> vsearch  --fastq_filter fastqfile  --fastq_maxee 1 --fastaout outputfile"),
                                                  hr(),
                                                  fluidRow(
                                                    column(
                                                      4,
                                                      h5(strong("Select the un-filtered reads folder (RST1_MergedFastq) to parse the directory...")),
                                                      shinyDirButton(
                                                        "rawUnfiltedReadsPath",
                                                        "Select Folder",
                                                        "Upload"
                                                      ),
                                                      h5(strong("Plot quality profiles of fastq-format reads and wait a few seconds 
                                                                (Support single file and =< 10 multiple files to display)...")),
                                                      shinyFilesButton(
                                                        "choseMergedFq",
                                                        "Select merged fastq file",
                                                        "Upload",
                                                        TRUE
                                                      ),
                                                      hr(),
                                                      actionButton('refreshQualityMerged',
                                                                   label = "Refresh Quality Profiles",
                                                                   icon = icon("refresh"))
                                                    ),
                                                    column(
                                                      5,
                                                      plotOutput(
                                                        outputId = "qualityMerged",
                                                        width = 400,height = 250,click = T
                                                      )
                                                      
                                                    )
                                                  ),
                                                  hr(),
                                                  h4("> Strip the specified number of bases from the left and right end of the reads"),
                                                  switchInput(
                                                    inputId = "stripLeftRightSwitch",
                                                    value = FALSE,
                                                    onStatus = "success",
                                                    offStatus = "danger",
                                                    size = "mini"
                                                  ),
                                                  fluidRow(
                                                    column(
                                                      4,
                                                      pickerInput(
                                                        inputId = "stripLeft",
                                                        label = "The number of bases from left end of read:",
                                                        choices = c(15, 16, 17, 18, 19, 20, 21, 22),
                                                        selected = 19
                                                      )
                                                    ),
                                                    column(
                                                      4,
                                                      pickerInput(
                                                        inputId = "stripRight",
                                                        label = "The number of bases from right end of read:",
                                                        choices = c(15, 16, 17, 18, 19, 20, 21, 22),
                                                        selected = 20
                                                      )
                                                    )
                                                  ),
                                                  hr(),
                                                  helpText("Next filte uploaded read files..."),
                                                  fluidRow(column(4, progressBar(id = "prgsBar_qualityFilter", value = 0, display_pct = T))),
                                                  actionBttn(
                                                    inputId = "qualityFilter",
                                                    label = "Quality filter reads...",
                                                    style = "jelly",
                                                    color = "primary"
                                                  )
                                                  ),
                                         tabPanel("DetectChimera",
                                                  helpText("* If  chimera reads (filted reads) were not detected, 
                                                           you should firstly upload that. 
                                                           File format: .fasta"),
                                                  helpText(">>Method-1: vsearch --uchime_ref fastafile  --db rdp_gold.fa --relabel Nochimera_ --nonchimeras outputfile"),
                                                  helpText(">>Method-2: vsearch --uchime_denovo fastafile --relabel Nochimera_ --nonchimeras outputfile"),
                                                  helpText(">>Method-3: vsearch --uchime2_denovo fastafile --relabel Nochimera_ --nonchimeras outputfile"),
                                                  helpText(">>Method-4: vsearch --uchime3_denovo fastafile --relabel Nochimera_ --nonchimeras outputfile"),
                                                  hr(),
                                                  fluidRow(
                                                    column(
                                                      3,
                                                      h5(strong("Select the undetected (filtered) reads folder (RST2_QualityFiltedFasta) ...")),
                                                      shinyDirButton(
                                                        "rawChimeraReadsPath",
                                                        "Select Folder",
                                                        "Upload"
                                                      )
                                                    ),
                                                    column(
                                                      3,
                                                      pickerInput(
                                                        inputId = "detectChimeraMethod",
                                                        label = "Select the method of detecting chimeras...",
                                                        choices = c("--uchime_ref", "--uchime_denovo", "--uchime2_denovo", "--uchime3_denovo"),
                                                        selected = "--uchime_ref"
                                                      )
                                                    )
                                                  ),
                                                  hr(),
                                                  helpText("Next detect chimera of uploaded read files..."),
                                                  fluidRow(column(4, progressBar(id = "prgsBar_chimeraDetect", value = 0, display_pct = T))),
                                                  actionBttn(
                                                    inputId = "chimeraDetect",
                                                    label = "Detect chimera...",
                                                    style = "jelly",
                                                    color = "primary"
                                                  )
                                                  ),
                                         tabPanel("Subsample",
                                                  helpText("* Each sample is reduced to a specified number of reads in order to reduce processing time. 
                                                           This subsampling allows you to ensure that the same number 
                                                           of reads are used from each sample in the OTU picking step. 
                                                           File format: .fasta"),
                                                  helpText(">> vsearch --fastx_subsample fastafile --randseed 123 --sample_size positive_integer --fastaout outputfile"),
                                                  hr(),
                                                  fluidRow(
                                                    column(
                                                      5,
                                                      h5(strong("Select the nochimeras reads folder (RST3_NoChimerasFasta)...")),
                                                      shinyDirButton(
                                                        "detectedChimeraReadsPath",
                                                        "Select Folder",
                                                        "Upload"
                                                      ),
                                                      h5(strong("Summary the number of sequences about all samples...")),
                                                      actionButton(
                                                        inputId = "summaryNumberSeqsDetectedChimeraFiles",
                                                        label = "Wait a few seconds", 
                                                        icon = icon("bar-chart-o")
                                                      ),
                                                      numericInput(
                                                        inputId = "setSubsampleSize",
                                                        label = "Set subsample size (Positive integer)...",
                                                        value = 200,
                                                        step = 1
                                                      )
                                                    ),
                                                    column(
                                                      5,
                                                      plotOutput(
                                                        outputId = "summaryDetectedChimeraSize",
                                                        width = 400,height = 250,click = T
                                                      )
                                                    )
                                                  ),
                                                  hr(),
                                                  helpText("Next subsample uploaded read files..."),
                                                  fluidRow(column(4, progressBar(id = "prgsBar_subsampleDetectedChimeraReads", value = 0, display_pct = T))),
                                                  actionBttn(
                                                    inputId = "subsampleDetectedChimeraReads",
                                                    label = "Subsample reads",
                                                    style = "jelly",
                                                    color = "primary"
                                                  )
                                                  ),
                                         tabPanel("Dereplicate",
                                                  helpText("* The dereplication step collapses all identical sequences to one 
                                                           and simply keeps track of how many there were. 
                                                           This can save a ton of time in further processing steps. 
                                                           If filted reads were not dereplicated, 
                                                           you should firstly upload that. 
                                                           File format: .fasta"),
                                                  helpText(">> vsearch --derep_fulllength fastafile --sizeout --minuniquesize positive_integer --output outputfile"),
                                                  hr(),
                                                  helpText("The nochimeras fastafiles or subsampled nochimeras fastafiles would be used..."),
                                                  fluidRow(
                                                    column(
                                                      3,
                                                      h5(strong("Upload undereplivcated reads folder (RST3_NoChimerasFasta or RST3_SubsampleNoChimerasFasta)...")),
                                                      shinyDirButton(
                                                        "rawUndereplivcatedReadsPath",
                                                        "Select Folder",
                                                        "Upload"
                                                      )
                                                    ),
                                                    column(
                                                      3,
                                                      pickerInput(
                                                        inputId = "selectMinuniquesize",
                                                        label = "Select min unique-size...",
                                                        choices = c(1,2,3,4,5),
                                                        selected = 2
                                                      )
                                                    )
                                                  ),
                                                  hr(),
                                                  helpText("Next dereplicate uploaded read files..."),
                                                  fluidRow(column(4, progressBar(id = "prgsBar_dereplicateReads", value = 0, display_pct = T))),
                                                  actionBttn(
                                                    inputId = "dereplicateReads",
                                                    label = "Dereplicate sequences...",
                                                    style = "jelly",
                                                    color = "primary"
                                                  )
                                                  )
                                         ) 
                                       ),
                                     tabPanel("Cluster OTU",
                                              br(),
                                              navbarPage(
                                                "Cluster sequences into OTUs >>",
                                                tabPanel("Cluster",
                                                         helpText("* If dereplicated reads were not clustered, you should firstly upload that. File format: .fasta "),
                                                         helpText(">> vsearch --cluster_fast fastafile --id realNum --iddef realNum --centroids OTU_fastafile --relabel OTU_"),
                                                         hr(),
                                                         fluidRow(
                                                           column(
                                                             4,
                                                             br(),br(),
                                                             h5(strong("Select dereplicated reads folder (RST4_DereplicatedFasta)...")),
                                                             shinyDirButton(
                                                               "rawUnclusterReadsPath",
                                                               "Select Folder",
                                                               "Upload"
                                                             )
                                                           ),
                                                           column(
                                                             4,
                                                             helpText("The pairwise identity is by default defined as 
                                                                      the number of (matching columns) / (alignment length - terminal gaps). 
                                                                      That definition can be modified by --iddef."),
                                                             sliderTextInput(
                                                               inputId = "selectClusterID",
                                                               label = "Select pairwise identity...",
                                                               choices = seq(0.7,1,0.01),
                                                               selected = 0.97,
                                                               grid = TRUE
                                                             )
                                                           )
                                                         ),
                                                         fluidRow(
                                                           column(
                                                             4,
                                                             h5("> Change the pairwise identity definition used in --id. Values accepted are:"),
                                                             h5("0. CD-HIT definition: (matching columns) / (shortest sequence length);"),
                                                             h5("1. edit distance: (matching columns) / (alignment length);"),
                                                             h5("2. edit distance excluding terminal gaps (same as --id);"),
                                                             h5("3. Marine Biological Lab definition counting each gap opening (internal orterminal) 
                                                                as a single mismatch, whether or not the gap was 
                                                                extended: 1.0 - [(mismatches + gap openings)/(longest sequence length)];"),
                                                             h5("4. BLAST definition, equivalent to --iddef 1 in a context of global pairwise alignment"),
                                                             br()
                                                           ),
                                                           column(
                                                             4,
                                                             br(),br(),br(),br(),br(),
                                                             radioGroupButtons(
                                                               inputId = "selectClusterIDdef",
                                                               label = "Select the pairwise identity definition (radio)...",
                                                               choices = c("ND",0,1,2,3,4),
                                                               checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                                                             )
                                                           )
                                                         ),
                                                         hr(),
                                                         helpText("Next cluster uploaded pair-reads..."),
                                                         fluidRow(column(4, progressBar(id = "prgsBar_clusterReads", value = 0, display_pct = T))),
                                                         actionBttn(
                                                           inputId = "clusterReads",
                                                           label = "Cluster dereplicated reads files...",
                                                           style = "jelly",
                                                           color = "primary"
                                                         )
                                                ),
                                                tabPanel("TrackReadsNumber",
                                                         helpText("Track number of reads at each step: 
                                                                  mergedfq, filtedfasta, deChimerafasta, subsample, dereplicated, OTUs. "),
                                                         actionButton(
                                                           inputId = "summaryPreparedReadsNumber",
                                                           label = "Wait a few seconds", 
                                                           icon = icon("bar-chart-o")
                                                         ),
                                                         hr(),
                                                         fluidRow(
                                                           column(5, progressBar(id = "prgsBar_summaryPreparedReadsNumber", value = 0, display_pct = T))
                                                         ),
                                                         hr(),
                                                         DT::dataTableOutput(outputId = "trackReadsNumber")
                                                )
                                              )
                                              ),
                                     tabPanel("OTU Table",
                                              br(),
                                              navbarPage(
                                                "Generate a count table >>",
                                                tabPanel("PerSample",
                                                         helpText("* The count table is what tells us how many times each sequence appears in each sample. 
                                                                  File format: .fasta "),
                                                         helpText(">> vsearch --usearch_global nochimerasFastafile --db OTUfastafile --id realNum --otutabout OTU_table"),
                                                         hr(),
                                                         fluidRow(
                                                           column(
                                                             4,
                                                             br(),
                                                             h5(strong("Select nochimeras reads folder (RST3_NoChimerasFasta or RST3_SubsampleNoChimerasFasta)...")),
                                                             shinyDirButton(
                                                               "rawnochimerasReadsPerPath",
                                                               "Select Folder",
                                                               "Upload"
                                                             )
                                                           ),
                                                           column(
                                                             4,
                                                             helpText("The pairwise identity is by default defined as 
                                                                      the number of (matching columns) / (alignment length - terminal gaps). 
                                                                      That definition can be modified by --iddef."),
                                                             sliderTextInput(
                                                               inputId = "selectConstructOTUidPer",
                                                               label = "Select pairwise identity...",
                                                               choices = seq(0.7,1,0.01),
                                                               selected = 0.97,
                                                               grid = TRUE
                                                             )
                                                           )
                                                         ),
                                                         fluidRow(
                                                           column(
                                                             4,
                                                             h5("> Change the pairwise identity definition used in --id. Values accepted are:"),
                                                             h5("0. CD-HIT definition: (matching columns) / (shortest sequence length);"),
                                                             h5("1. edit distance: (matching columns) / (alignment length);"),
                                                             h5("2. edit distance excluding terminal gaps (same as --id);"),
                                                             h5("3. Marine Biological Lab definition counting each gap opening (internal orterminal) 
                                                                as a single mismatch, whether or not the gap was 
                                                                extended: 1.0 - [(mismatches + gap openings)/(longest sequence length)];"),
                                                             h5("4. BLAST definition, equivalent to --iddef 1 in a context of global pairwise alignment"),
                                                             br()
                                                             ),
                                                           column(
                                                             4,
                                                             br(),br(),br(),br(),
                                                             radioGroupButtons(
                                                               inputId = "selectOTUsTablePerIDdef",
                                                               label = "Select the pairwise identity definition (radio)...",
                                                               choices = c("ND",0,1,2,3,4),
                                                               checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                                                             )
                                                           )
                                                           ),
                                                         hr(),
                                                         helpText("Next summary OTU table per sample..."),
                                                         fluidRow(column(4, progressBar(id = "prgsBar_summaryOTUtablePer", value = 0, display_pct = T))),
                                                         actionBttn(
                                                           inputId = "summaryOTUtablePer",
                                                           label = "Summary OTU table per sample...",
                                                           style = "jelly",
                                                           color = "primary"
                                                         )
                                                         ),
                                                tabPanel("AllSamples",
                                                         helpText("* The count table is what tells us how many 
                                                                  times each sequence appears among all samples. File format: .fasta"),
                                                         hr(),
                                                         fluidRow(
                                                           column(
                                                             4,
                                                             br(),
                                                             h5(strong("Select nochimeras reads folder (RST3_NoChimerasFasta or RST3_SubsampleNoChimerasFasta) ...")),
                                                             shinyDirButton(
                                                               "rawnochimerasReadsAllPath",
                                                               "Select Folder",
                                                               "Upload"
                                                             )
                                                           ),
                                                           column(
                                                             4,
                                                             helpText("The pairwise identity is by default defined as 
                                                                      the number of (matching columns) / (alignment length - terminal gaps). 
                                                                      That definition can be modified by --iddef."),
                                                             sliderTextInput(
                                                               inputId = "selectConstructOTUidAll",
                                                               label = "Select pairwise identity...",
                                                               choices = seq(0.7,1,0.01),
                                                               selected = 0.97,
                                                               grid = TRUE
                                                             )
                                                           )
                                                         ),
                                                         fluidRow(
                                                           column(
                                                             4,
                                                             h5("> Change the pairwise identity definition used in --id. Values accepted are:"),
                                                             h5("0. CD-HIT definition: (matching columns) / (shortest sequence length);"),
                                                             h5("1. edit distance: (matching columns) / (alignment length);"),
                                                             h5("2. edit distance excluding terminal gaps (same as --id);"),
                                                             h5("3. Marine Biological Lab definition counting each gap opening (internal orterminal) 
                                                                as a single mismatch, whether or not the gap was 
                                                                extended: 1.0 - [(mismatches + gap openings)/(longest sequence length)];"),
                                                             h5("4. BLAST definition, equivalent to --iddef 1 in a context of global pairwise alignment"),
                                                             br()
                                                           ),
                                                           column(
                                                             4,
                                                             br(),br(),br(),
                                                             radioGroupButtons(
                                                               inputId = "selectOTUsTableAllIDdef",
                                                               label = "Select the pairwise identity definition (radio)...",
                                                               choices = c("ND",0,1,2,3,4),
                                                               checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                                                             ),
                                                             br()
                                                           )
                                                         ),
                                                         hr(),
                                                         helpText("Next summary OTU table of all samples..."),
                                                         fluidRow(
                                                           column(6, 
                                                                  progressBar(id = "prgsBar_summaryNochimerasFastaAll", 
                                                                              title = "Process-1",
                                                                              value = 0, 
                                                                              display_pct = T),
                                                                  progressBar(id = "prgsBar_summaryTotalOTUsFastaAll", 
                                                                              title = "Process-2",
                                                                              value = 0, 
                                                                              display_pct = T),
                                                                  progressBar(id = "prgsBar_summaryAllOTUtableEach", 
                                                                              title = "Process-3",
                                                                              value = 0, 
                                                                              display_pct = T),
                                                                  progressBar(id = "prgsBar_summaryAllOTUsMatrix", 
                                                                              title = "Process-4",
                                                                              value = 0, 
                                                                              display_pct = T))
                                                           ),
                                                         actionBttn(
                                                           inputId = "summaryOTUtableAll",
                                                           label = "Summary OTU table of all samples...",
                                                           style = "jelly",
                                                           color = "primary"
                                                         ),
                                                         hr(),
                                                         DT::dataTableOutput(outputId = "showAllOTUsMatrix")
                                                         )
                                              )
                                              ),
                                     tabPanel("Assign Taxonomy",
                                              br(),
                                              navbarPage(
                                                "Taxonomic classification >>",
                                                tabPanel("Method-1",
                                                         helpText("* The final step in our sequence processing is to assign taxonomy to OTUs. "),
                                                         helpText(">> vsearch --usearch_global OTUfastafile --db database --id realNum --blast6out taxOutput "),
                                                         hr(),
                                                         fluidRow(
                                                           column(
                                                             4,
                                                             br(),
                                                             materialSwitch(
                                                               inputId = "assignMethod1Status",
                                                               label = "Use Method-1 ? (Default: OFF) ", 
                                                               value = FALSE,
                                                               status = "success"
                                                             ),
                                                             br(),
                                                             checkboxGroupButtons(
                                                               inputId = "selectDBstyleGlobal",
                                                               label = "Select the type of ribosomal database...",
                                                               choices = c("SILVA", "RDP", "GreenGene", "EzBioCloud"),
                                                               selected = c("SILVA"),
                                                               checkIcon = list(
                                                                 yes = tags$i(class = "fa fa-check-square", 
                                                                              style = "color: steelblue"),
                                                                 no = tags$i(class = "fa fa-square-o", 
                                                                             style = "color: steelblue"))
                                                             )
                                                             ),
                                                           column(
                                                             4,
                                                             helpText("The pairwise identity is by default defined as 
                                                                      the number of (matching columns) / (alignment length - terminal gaps). 
                                                                      That definition can be modified by --iddef."),
                                                             sliderTextInput(
                                                               inputId = "selectGlobalID",
                                                               label = "Select pairwise identity...",
                                                               choices = seq(0.7,1,0.01),
                                                               selected = 0.97,
                                                               grid = TRUE
                                                             )
                                                             )
                                                         ),
                                                         fluidRow(
                                                           column(
                                                             4,
                                                             h5("> Change the pairwise identity definition used in --id. Values accepted are:"),
                                                             h5("0. CD-HIT definition: (matching columns) / (shortest sequence length);"),
                                                             h5("1. edit distance: (matching columns) / (alignment length);"),
                                                             h5("2. edit distance excluding terminal gaps (same as --id);"),
                                                             h5("3. Marine Biological Lab definition counting each gap opening (internal orterminal) 
                                                                as a single mismatch, whether or not the gap was 
                                                                extended: 1.0 - [(mismatches + gap openings)/(longest sequence length)];"),
                                                             h5("4. BLAST definition, equivalent to --iddef 1 in a context of global pairwise alignment")
                                                           ),
                                                           column(
                                                             4,
                                                             br(),br(),br(),br(),br(),
                                                             checkboxGroupButtons(
                                                               inputId = "selectGlobalIDdef",
                                                               label = "Select the pairwise identity definition (radio)...",
                                                               choices = c(0,1,2,3,4),
                                                               checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                                                             )
                                                           )
                                                         ),
                                                         br(),hr(),
                                                         fluidRow(
                                                           column(
                                                             4,
                                                             h4("> Step-1: assign taxonomy to OTUs of each sample"),
                                                             fileInput(
                                                               inputId = "rawOTUfastaPerGlobalPath",
                                                               label = "Select anyone local OTUs fasta file of each sample to parse the directory...",
                                                               buttonLabel = "Upload...",
                                                               multiple = FALSE
                                                             ),
                                                             helpText("Next assign taxonomy uploaded OTUs fasta..."),
                                                             fluidRow(column(8, progressBar(id = "prgsBar_assignTax2OTUperGlobal", value = 0, display_pct = T))),
                                                             actionBttn(
                                                               inputId = "assignTax2OTUperGlobal",
                                                               label = "Assign taxonomy...",
                                                               style = "jelly",
                                                               color = "primary"
                                                             )
                                                           ),
                                                           column(
                                                             4,
                                                             h4("> Step-2: assign taxonomy to summaried-OTUs of all samples"),
                                                             fileInput(
                                                               inputId = "rawOTUfastaAllGlobalPath",
                                                               label = "Select anyone local OTUs fasta files summaried from all samples to parse the directory...",
                                                               buttonLabel = "Upload...",
                                                               multiple = FALSE
                                                             ),
                                                             helpText("Next assign taxonomy uploaded summaried-OTUs fasta..."),
                                                             fluidRow(column(8, progressBar(id = "prgsBar_assignTax2OTUallGlobal", value = 0, display_pct = T))),
                                                             actionBttn(
                                                               inputId = "assignTax2OTUallGlobal",
                                                               label = "Assign taxonomy...",
                                                               style = "jelly",
                                                               color = "primary"
                                                             )
                                                           )
                                                         )
                                                ),
                                                tabPanel("Method-2",
                                                         helpText("* The final step in our sequence processing is to assign taxonomy to OTUs. "),
                                                         helpText(">> vsearch --sintax OTUfastafile --db database --sintax_cutoff realNum --tabbedout taxOutput"),
                                                         hr(),
                                                         fluidRow(
                                                           column(
                                                             4,
                                                             br(),
                                                             materialSwitch(
                                                               inputId = "assignMethod2Status",
                                                               label = "Use Method-2 ? (Default: OFF) ", 
                                                               value = FALSE,
                                                               status = "success"
                                                             ),
                                                             br(),
                                                             checkboxGroupButtons(
                                                               inputId = "selectDBstyleSintax",
                                                               label = "Select the type of ribosomal database...",
                                                               choices = c("SILVA", "RDP", "GreenGene", "EzBioCloud"),
                                                               selected = c("SILVA"),
                                                               checkIcon = list(
                                                                 yes = tags$i(class = "fa fa-check-square", 
                                                                              style = "color: steelblue"),
                                                                 no = tags$i(class = "fa fa-square-o", 
                                                                             style = "color: steelblue"))
                                                             )
                                                           ),
                                                           column(
                                                             4,
                                                             helpText("The --sintax_cutoff option may be used to set a minimum 
                                                                      level of bootstrap support for the taxonomic ranks to be reported."),
                                                             sliderTextInput(
                                                               inputId = "selectSintaxCutoff",
                                                               label = "Select sintax cutoff...",
                                                               choices = seq(0.7,1,0.01),
                                                               selected = 0.97,
                                                               grid = TRUE
                                                             )
                                                             )
                                                         ),
                                                         hr(),
                                                         fluidRow(
                                                           column(
                                                             4,
                                                             h4("> Step-1: assign taxonomy to OTUs of each sample"),
                                                             fileInput(
                                                               inputId = "rawOTUfastaPerSintaxPath",
                                                               label = "Select anyone local OTUs fasta file of each sample to parse the directory...",
                                                               buttonLabel = "Upload...",
                                                               multiple = FALSE
                                                             ),
                                                             helpText("Next assign taxonomy uploaded OTUs fasta..."),
                                                             fluidRow(column(8, progressBar(id = "prgsBar_assignTax2OTUperSintax", value = 0, display_pct = T))),
                                                             actionBttn(
                                                               inputId = "assignTax2OTUperSintax",
                                                               label = "Assign taxonomy...",
                                                               style = "jelly",
                                                               color = "primary"
                                                             )
                                                           ),
                                                           column(
                                                             4,
                                                             h4("> Step-2: assign taxonomy to summaried-OTUs of all samples"),
                                                             fileInput(
                                                               inputId = "rawOTUfastaAllSintaxPath",
                                                               label = "Select anyone local OTUs fasta files summaried from all samples to parse the directory...",
                                                               buttonLabel = "Upload...",
                                                               multiple = FALSE
                                                             ),
                                                             helpText("Next assign taxonomy uploaded summaried-OTUs fasta..."),
                                                             fluidRow(column(8, progressBar(id = "prgsBar_assignTax2OTUallSintax", value = 0, display_pct = T))),
                                                             actionBttn(
                                                               inputId = "assignTax2OTUallSintax",
                                                               label = "Assign taxonomy...",
                                                               style = "jelly",
                                                               color = "primary"
                                                             )
                                                           )
                                                         )
                                                         )
                                              )
                                              ),
                                     tabPanel("Analysis Diversity",
                                              br(),
                                              navbarPage(
                                                "Analysising diversities >>",
                                                tabPanel("α-diversity",
                                                         helpText("
                                                           * If pair-reads were not merged, 
                                                           you should firstly upload that. 
                                                           File format: .fq or .fastq. "),
                                                         br(),
                                                         fileInput(
                                                           inputId = "rawPairReads",
                                                           label = "Select local pair-reads files...",
                                                           buttonLabel = "Upload...",
                                                           multiple = TRUE
                                                         ),
                                                         hr(),
                                                         helpText("Next merge uploaded pair-reads..."),
                                                         fluidRow(column(4, progressBar(id = "prgsBar_mergePairs", value = 0, display_pct = T))),
                                                         actionBttn(
                                                           inputId = "mergePairs",
                                                           label = "Merge Pair-reads files...",
                                                           style = "jelly",
                                                           color = "primary"
                                                         )
                                                         
                                                ),
                                                tabPanel("β-diversity",
                                                         helpText("* If pair-reads were not merged, you should firstly upload that. File format: .fq or .fastq. "),
                                                         br(),
                                                         fileInput(
                                                           inputId = "rawPairReads",
                                                           label = "Select local pair-reads files...",
                                                           buttonLabel = "Upload...",
                                                           multiple = TRUE
                                                         ),
                                                         hr(),
                                                         helpText("Next merge uploaded pair-reads..."),
                                                         fluidRow(column(4, progressBar(id = "prgsBar_mergePairs", value = 0, display_pct = T))),
                                                         actionBttn(
                                                           inputId = "mergePairs",
                                                           label = "Merge Pair-reads files...",
                                                           style = "jelly",
                                                           color = "primary"
                                                         )
                                                )
                                              )
                                              )
                                     )
                                   )
                                 ),
                        
                        
                        
                        tabPanel("BGDstudio")
                        
                        
                        
                        
                        
                        )
            )



server <- function(input, output, session) {
  
  # extract full path of selected folder
  getFullPath <- function(selectFolder){
    dirUnlsit <- as.vector(unlist(selectFolder))
    fullpath <- c()
    rootpath <- strsplit(strsplit(dirUnlsit[length(dirUnlsit)],'(',fixed = T)[[1]],')',fixed = T)[[2]]
    fullpath <- paste0(rootpath, '/')
    for (pathNum in seq(2:length(dirUnlsit))) {
      if (length(dirUnlsit) - pathNum >= 1 & dirUnlsit[pathNum] > 1 ) {
        fullpath <- paste0(fullpath,dirUnlsit[pathNum], '/')
      }
    }
    return(fullpath)
  }
  
  # extract parent full path of selected folder
  getParentFullPath <- function(selectFolder){
    dirUnlsit <- as.vector(unlist(selectFolder))
    fullpath <- c()
    rootpath <- strsplit(strsplit(dirUnlsit[length(dirUnlsit)],'(',fixed = T)[[1]],')',fixed = T)[[2]]
    fullpath <- paste0(rootpath, '/')
    for (pathNum in seq(2:length(dirUnlsit))) {
      if (length(dirUnlsit) - pathNum >= 2 & dirUnlsit[pathNum] > 1 ) {
        fullpath <- paste0(fullpath,dirUnlsit[pathNum], '/')
      }
    }
    return(fullpath)
  }
  
  # extract OTUs table from 'construct OTUs Table for each sample' & rewrite the --usearch_global outoupt 
  extractOTUsTable <- function(folderPath, inputFilename){
    file.df <- read.table( str_c(folderPath, inputFilename) , 
                           header = F,sep = '\t',row.names = 1, 
                           stringsAsFactors = F)
    OTUs.DF <- data.frame(OTUs = rownames(file.df), Counter = rowSums(file.df))
    
    write.table(OTUs.DF, str_c(folderPath, inputFilename), 
                quote = F,sep = '\t',row.names = F, col.names = F)
  }
  
  volumes = getVolumes()()
  
  
  #######################################################################################
  ####    Prepare Reads    ##############################################################
  #######################################################################################
  # ==================================================================
  # Merge Pair-end reads
  shinyDirChoose(input,
                 'rawPairReadsPath',
                 roots = volumes)
  
  observeEvent(input$mergePairs, {
    parentPEfilesPath <- getParentFullPath(input$rawPairReadsPath)
    rawPEfilesPath <- getFullPath(input$rawPairReadsPath)
    symbolForward <- as.character(input$difsymbolForward)
    symbolReverse <- as.character(input$difsymbolReverse)
    sampleNames <- c()
    for (filename in list.files(rawPEfilesPath)) {
      if ( str_detect(filename, symbolForward) ) {
        sampleNames <- append(sampleNames, strsplit(filename, symbolForward)[[1]])
      }
    }
    
    # creat & refresh RST1_MergedFastq
    if (dir.exists(str_c(parentPEfilesPath,'RST1_MergedFastq/'))) {
      unlink(str_c(parentPEfilesPath,'RST1_MergedFastq/'), recursive = TRUE)
    } else {
      dir.create(str_c(parentPEfilesPath,'RST1_MergedFastq/'))
    }
    
    # start to merge pair-end files
    mergeComd1 <- 'vsearch.exe --fastq_mergepairs '
    mergeComd2 <- ' --reverse '
    mergeComd3 <- ' --fastqout '
    sampleNum <- length(sampleNames)
    for (i in 1:sampleNum) {
      samplename <- sampleNames[i]
      forwardfq <- str_c(rawPEfilesPath, samplename, symbolForward)
      reversefq <- str_c(rawPEfilesPath, samplename, symbolReverse)
      fastqout <- str_c(parentPEfilesPath, 'RST1_MergedFastq/', samplename, '.fastq')
      system(str_c(mergeComd1, forwardfq, mergeComd2, reversefq, mergeComd3, fastqout))
      updateProgressBar(session, id = "prgsBar_mergePairs", value = (i/sampleNum)*100)
    }
  })
  
  
  # ==================================================================
  # Quality-filter
  shinyDirChoose(input,
                 'rawUnfiltedReadsPath',
                 roots = volumes)
  
  observe({
    shinyFileChoose(input,
                    "choseMergedFq",
                    roots = volumes,
                    session = session)
    observeEvent(input$choseMergedFq,{
      if(!is.null(input$choseMergedFq)){
        file_selected<-parseFilePaths(volumes, input$choseMergedFq)
        output$qualityMerged <- renderPlot({
          plotQualityProfile(as.character(file_selected$datapath)[1])
        })
      }
    })
  
    refreshClick <<- 1
    observeEvent(input$refreshQualityMerged,{
      file_selected<-parseFilePaths(volumes, input$choseMergedFq)
      length(file_selected$datapath)
      if (length(file_selected$datapath) >= 2) {
        output$qualityMerged <- renderPlot({
          plotQualityProfile(as.character(file_selected$datapath)[refreshClick])
        })
      }
      refreshClick <<- refreshClick + 1
      if (length(file_selected$datapath) < refreshClick ){
        output$qualityMerged <- renderPlot({
          barplot(1,legend.text = NULL,
                  ylab = NULL,xlab = NULL,
                  axes = F,col = "white", col.main = 'red',
                  main = "Out Range of Selected Files")})
      }
    })
  })
  
  observeEvent(input$qualityFilter, {
    parentMergedfilesPath <- getParentFullPath(input$rawUnfiltedReadsPath)
    rawMergedfilesPath <- getFullPath(input$rawUnfiltedReadsPath)
    sampleNames <- c()
    for (filename in list.files(rawMergedfilesPath)) {
      sampleNames <- append(sampleNames, strsplit(filename, '.fastq')[[1]])
    }
    
    # creat & refresh RST2_QualityFiltedFasta
    if (dir.exists(str_c(parentMergedfilesPath,'RST2_QualityFiltedFasta/'))) {
      unlink(str_c(parentMergedfilesPath,'RST2_QualityFiltedFasta/'), recursive = TRUE)
    } else {
      dir.create(str_c(parentMergedfilesPath,'RST2_QualityFiltedFasta/'))
    }
    
    # start to filter merged files
    filtedComd1 <- 'vsearch.exe --fastq_filter '
    filtedComd2 <- ' --fastq_maxee 1 --fastaout '
    sampleNum <- length(sampleNames)
    if (input$stripLeftRightSwitch == T) {
      input$stripLeft
      input$stripRight
      filtedComd3 <- str_c(' --fastq_stripleft ', input$stripLeft, ' --fastq_stripright ', input$stripRight)
      for (i in 1:sampleNum) {
        samplename <- sampleNames[i]
        mergedfq <- str_c(rawMergedfilesPath, samplename, '.fastq')
        fastaout <- str_c(parentMergedfilesPath, 'RST2_QualityFiltedFasta/', samplename, '.fasta')
        system(str_c(filtedComd1, mergedfq, filtedComd2,  fastaout, filtedComd3))
        updateProgressBar(session, id = "prgsBar_qualityFilter", value = (i/sampleNum)*100)
      }
    } else {
      for (i in 1:sampleNum) {
        samplename <- sampleNames[i]
        mergedfq <- str_c(rawMergedfilesPath, samplename, '.fastq')
        fastaout <- str_c(parentMergedfilesPath, 'RST2_QualityFiltedFasta/', samplename, '.fasta')
        system(str_c(filtedComd1, mergedfq, filtedComd2,  fastaout))
        updateProgressBar(session, id = "prgsBar_qualityFilter", value = (i/sampleNum)*100)
      }
    }
  })
  
  
  # ==================================================================
  # Detect Chimeras
  shinyDirChoose(input,
                 'rawChimeraReadsPath',
                 roots = volumes)
  
  observeEvent(input$chimeraDetect, {
    parentFilteredfilesPath <- getParentFullPath(input$rawChimeraReadsPath)
    rawChimerafilesPath <- getFullPath(input$rawChimeraReadsPath)
    sampleNames <- list.files(rawChimerafilesPath)
    
    # creat & refresh RST3_NoChimerasFasta
    if (dir.exists(str_c(parentFilteredfilesPath,'RST3_NoChimerasFasta/'))) {
      unlink(str_c(parentFilteredfilesPath,'RST3_NoChimerasFasta/'), recursive = TRUE)
    } else {
      dir.create(str_c(parentFilteredfilesPath,'RST3_NoChimerasFasta/'))
    }
    
    # start to detect chimeras of filtered files
    filtedComd1 <- str_c('vsearch.exe  --relabel Nochimera_ --nonchimeras ', parentFilteredfilesPath,'RST3_NoChimerasFasta/')
    sampleNum <- length(sampleNames)
    
    if (input$detectChimeraMethod == '--uchime_ref') {
      filtedComd2 <- str_c(' --db ', 'rdp_gold.fasta ')
      filtedComd3 <- str_c(' --uchime_ref ', rawChimerafilesPath)
      for (i in 1:sampleNum) {
        samplename <- sampleNames[i]
        system(str_c(filtedComd1, samplename, filtedComd2, filtedComd3, samplename))
        updateProgressBar(session, id = "prgsBar_chimeraDetect", value = (i/sampleNum)*100)
      }
    } else if (input$detectChimeraMethod == '--uchime_denovo') {
      filtedComd2 <- str_c(' --uchime_denovo ', rawChimerafilesPath)
      for (i in 1:sampleNum) {
        samplename <- sampleNames[i]
        system(str_c(filtedComd1, samplename, filtedComd2, samplename))
        updateProgressBar(session, id = "prgsBar_chimeraDetect", value = (i/sampleNum)*100)
      }
    } else if (input$detectChimeraMethod == '--uchime2_denovo') {
      filtedComd2 <- str_c(' --uchime2_denovo ', rawChimerafilesPath)
      for (i in 1:sampleNum) {
        samplename <- sampleNames[i]
        system(str_c(filtedComd1, samplename, filtedComd2, samplename))
        updateProgressBar(session, id = "prgsBar_chimeraDetect", value = (i/sampleNum)*100)
      }
    } else if (input$detectChimeraMethod == '--uchime3_denovo') {
      filtedComd2 <- str_c(' --uchime3_denovo ', rawChimerafilesPath)
      for (i in 1:sampleNum) {
        samplename <- sampleNames[i]
        system(str_c(filtedComd1, samplename, filtedComd2, samplename))
        updateProgressBar(session, id = "prgsBar_chimeraDetect", value = (i/sampleNum)*100)
      }
    }
  })
  
  
  # ==================================================================
  # subsample
  shinyDirChoose(input,
                 'detectedChimeraReadsPath',
                 roots = volumes)
  
  # plot
  observeEvent(input$summaryNumberSeqsDetectedChimeraFiles, {
    output$summaryDetectedChimeraSize <- renderPlot({
      parentNochimerasfilesPath <- getParentFullPath(input$detectedChimeraReadsPath)
      rawNochimerasfilesPath <- getFullPath(input$detectedChimeraReadsPath)
      sampleNames <- c()
      seqsSummary <- c()
      for (filename in list.files(rawNochimerasfilesPath)) {
        sampleNames <- append(sampleNames, strsplit(filename, '.fasta')[[1]])
        seqsNum <- as.integer(length(getSequences(str_c(rawNochimerasfilesPath, filename))))
        seqsSummary <- append(seqsSummary, seqsNum)
      }
      seqSummaryDf <- data.frame(samples = sampleNames, seqsNum = seqsSummary)
      ggplot(seqSummaryDf, aes(x=samples, y=seqsNum)) +
        geom_point(size=2) +
        geom_segment(aes(x=samples, xend=samples, y=0, yend=seqsNum)) +
        geom_hline(aes(yintercept = min(seqsSummary)), colour = "red", linetype="dashed") +
        labs(title="Summary the number of sequences among all samples",
             subtitle= str_c('The min size - - ', min(seqsSummary),
                             ' ; Total samples - - ', length(seqsSummary))) +
        xlab(label = "Total samples") + ylab(label = "The number of sequences") +
        theme_bw() +theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1))
    })
  })
  
  # start to subsample
  observeEvent(input$subsampleDetectedChimeraReads, {
    parentNochimerasfilesPath <- getParentFullPath(input$detectedChimeraReadsPath)
    rawNochimerasfilesPath <- getFullPath(input$detectedChimeraReadsPath)
    sampleNames <- list.files(rawNochimerasfilesPath)

    # creat & refresh RST3_SubsampleNoChimerasFasta
    if (dir.exists(str_c(parentNochimerasfilesPath,'RST3_SubsampleNoChimerasFasta/'))) {
      unlink(str_c(parentNochimerasfilesPath,'RST3_SubsampleNoChimerasFasta/'), recursive = TRUE)
    } else {
      dir.create(str_c(parentNochimerasfilesPath,'RST3_SubsampleNoChimerasFasta/'))
    }

    # start to subsample nochimeras files
    subsampleComd1 <- str_c('python ./Scripts/Subsample.py ',input$setSubsampleSize, ' ')
    subsampleComd2 <- str_c(' ', parentNochimerasfilesPath, 'RST3_SubsampleNoChimerasFasta/')
    sampleNum <- length(sampleNames)
    for (i in 1:sampleNum) {
      samplename <- sampleNames[i]
      system(str_c(subsampleComd1, rawNochimerasfilesPath, samplename, subsampleComd2, samplename))
      updateProgressBar(session, id = "prgsBar_subsampleDetectedChimeraReads", value = (i/sampleNum)*100)
    }
  })
  
  
  # ==================================================================
  # Dereplication
  shinyDirChoose(input,
                 'rawUndereplivcatedReadsPath',
                 roots = volumes)
  
  observeEvent(input$dereplicateReads, {
    parentnochimerasfilesPath <- getParentFullPath(input$rawUndereplivcatedReadsPath)
    rawnochimerasfilesPath <- getFullPath(input$rawUndereplivcatedReadsPath)
    sampleNames <- list.files(rawnochimerasfilesPath)
    sampleNum <- length(sampleNames)
    
    # creat & refresh RST4_DereplicatedFasta
    if (dir.exists(str_c(parentnochimerasfilesPath,'RST4_DereplicatedFasta/'))) {
      unlink(str_c(parentnochimerasfilesPath,'RST4_DereplicatedFasta/'), recursive = TRUE)
    } else {
      dir.create(str_c(parentnochimerasfilesPath,'RST4_DereplicatedFasta/'))
    }
    
    # start to dereplicate nochimeras files
    dereplicateComd1 <- str_c('vsearch.exe --derep_fulllength ', rawnochimerasfilesPath)
    dereplicateComd2 <- str_c(' --sizeout --minuniquesize ', input$selectMinuniquesize)
    dereplicateComd3 <- str_c(' --output ', parentnochimerasfilesPath, 'RST4_DereplicatedFasta/')
    for (i in 1:sampleNum) {
      samplename <- sampleNames[i]
      system(str_c(dereplicateComd1, samplename, dereplicateComd2, dereplicateComd3, samplename))
      updateProgressBar(session, id = "prgsBar_dereplicateReads", value = (i/sampleNum)*100)
    }
  })
  
  
  #######################################################################################
  ####    Cluster OTU    ################################################################
  #######################################################################################
  # ==================================================================
  # cluster OTU
  shinyDirChoose(input,
                 'rawUnclusterReadsPath',
                 roots = volumes)
  
  observeEvent(input$clusterReads, {
    parentNoclusterfilesPath <- getParentFullPath(input$rawUnclusterReadsPath)
    rawNoclusterfilesPath <- getFullPath(input$rawUnclusterReadsPath)
    sampleNames <- list.files(rawNoclusterfilesPath)
    sampleNum <- length(sampleNames)
    
    # creat & refresh RST5_ClusterOTUfasta
    if (dir.exists(str_c(parentNoclusterfilesPath,'RST5_ClusterOTUfasta/'))) {
      unlink(str_c(parentNoclusterfilesPath,'RST5_ClusterOTUfasta/'), recursive = TRUE)
    } else {
      dir.create(str_c(parentNoclusterfilesPath,'RST5_ClusterOTUfasta/'))
    }
    
    # start to cluster dereplicated files
    clusterComd1 <- str_c('vsearch.exe --cluster_fast ', rawNoclusterfilesPath)
    clusterComd2 <- str_c(' --id ', input$selectClusterID)
    clusterComd3 <- str_c(' --relabel OTU_  --centroids ', parentNoclusterfilesPath, 'RST5_ClusterOTUfasta/' )
    clusterComd4 <- str_c(' --iddef ', input$selectClusterIDdef)
    for (i in 1:sampleNum) {
      samplename <- sampleNames[i]
      if (input$selectClusterIDdef == "ND") {
        system(str_c(clusterComd1, samplename, clusterComd2, clusterComd3, samplename))
        updateProgressBar(session, id = "prgsBar_clusterReads", value = (i/sampleNum)*100)
      } else {
        system(str_c(clusterComd1, samplename, clusterComd2, clusterComd3, samplename, clusterComd4))
        updateProgressBar(session, id = "prgsBar_clusterReads", value = (i/sampleNum)*100)
      }
    }
  })
  
  observeEvent(input$summaryPreparedReadsNumber, {
    parentNoclusterfilesPath <- getParentFullPath(input$rawUnclusterReadsPath)
    rawNoclusterfilesPath <- getFullPath(input$rawUnclusterReadsPath)
    setwd(parentNoclusterfilesPath)
    sampleNames <- c()
    for (fn in list.files(rawNoclusterfilesPath)) {
      sampleNames <- append(sampleNames, strsplit(fn, '.fasta')[[1]] )
    }
    # summary the number of mergedFq, filteredFasta, noChimeraFasta, subsampleFasta, dereplicatedFasta, and OTUsFasta
    getReadsNumVector <- function(seqFilesPath){
      if (file.exists(str_c('./', seqFilesPath)) & length(list.files(str_c('./', seqFilesPath))) != 0) {
        seqFilesPath.files <- list.files(seqFilesPath)
        seqFilesPath.readsNum <- c()
        for (fn in seqFilesPath.files) {
          seqsNum <- as.integer(length(getSequences(str_c(seqFilesPath, '/', fn))))
          seqFilesPath.readsNum <- append(seqFilesPath.readsNum, seqsNum)
        }
      } else {
        seqFilesPath.readsNum <- c()
      }
      return(seqFilesPath.readsNum)
    }
    
    updateProgressBar(session, id = "prgsBar_summaryPreparedReadsNumber", value = 5)
    mergedFq.readsNum <- getReadsNumVector('RST1_MergedFastq')
    updateProgressBar(session, id = "prgsBar_summaryPreparedReadsNumber", value = 20)
    filteredFasta.readsNum <- getReadsNumVector('RST2_QualityFiltedFasta')
    updateProgressBar(session, id = "prgsBar_summaryPreparedReadsNumber", value = 40)
    noChimeraFasta.readsNum <- getReadsNumVector('RST3_NoChimerasFasta')
    updateProgressBar(session, id = "prgsBar_summaryPreparedReadsNumber", value = 60)
    subsampleFasta.readsNum <- getReadsNumVector('RST3_SubsampleNoChimerasFasta')
    updateProgressBar(session, id = "prgsBar_summaryPreparedReadsNumber", value = 80)
    dereplicatedFasta.readsNum <- getReadsNumVector('RST4_DereplicatedFasta')
    updateProgressBar(session, id = "prgsBar_summaryPreparedReadsNumber", value = 90)
    OTUsFasta.readsNum <- getReadsNumVector('RST5_ClusterOTUfasta')
    updateProgressBar(session, id = "prgsBar_summaryPreparedReadsNumber", value = 100)
    output$trackReadsNumber <- DT::renderDataTable({
      trackDF <- data.frame(mergedFq.readsNum, filteredFasta.readsNum, 
                            noChimeraFasta.readsNum, subsampleFasta.readsNum, 
                            dereplicatedFasta.readsNum, OTUsFasta.readsNum)
      colnames(trackDF) <- c('mergedFq', 'filteredFast', 'noChimeraFasta', 
                             'subsampleFasta', 'dereplicatedFast', 'OTUsFasta')
      rownames(trackDF) <- sampleNames
      trackDF
    })
    
  })
  
  
  # ==================================================================
  # Construct OTUs Table -- for each sample
  shinyDirChoose(input,
                 'rawnochimerasReadsPerPath',
                 roots = volumes)
  
  observeEvent(input$summaryOTUtablePer, {
    parentnochimerasPathper <- getParentFullPath(input$rawnochimerasReadsPerPath)
    rawnochimerassPathper <- getFullPath(input$rawnochimerasReadsPerPath)
    sampleNames <- list.files(rawnochimerassPathper)
    sampleNum <- length(sampleNames)
    
    # creat & refresh RST6_OTUsTablePerSample
    if (dir.exists(str_c(parentnochimerasPathper,'RST6_OTUsTablePerSample/'))) {
      unlink(str_c(parentnochimerasPathper,'RST6_OTUsTablePerSample/'), recursive = TRUE)
    } else {
      dir.create(str_c(parentnochimerasPathper,'RST6_OTUsTablePerSample/'))
    }
    
    # start to construct OTUs table for each sample
    perTableOTUsComd1 <- str_c('vsearch.exe --usearch_global ', rawnochimerassPathper)
    perTableOTUsComd2 <- str_c(' --id ', input$selectConstructOTUidPer)
    perTableOTUsComd3 <- str_c(' --db ', parentnochimerasPathper,'RST5_ClusterOTUfasta/')
    perTableOTUsComd4 <- str_c(' --otutabout ', parentnochimerasPathper, 'RST6_OTUsTablePerSample/')
    perTableOTUsComd5 <- str_c(' --iddef ', input$selectOTUsTablePerIDdef)
    for (i in 1:sampleNum) {
      samplename <- sampleNames[i]
      otuoutputFileName <- str_c(strsplit(samplename, '.fasta')[[1]], '.txt')
      if (input$selectOTUsTablePerIDdef == "ND") {
        system(str_c(perTableOTUsComd1, samplename, perTableOTUsComd2, perTableOTUsComd3, 
                     samplename, perTableOTUsComd4, otuoutputFileName))
        extractOTUsTable(str_c(parentnochimerasPathper, 'RST6_OTUsTablePerSample/'), otuoutputFileName)
        updateProgressBar(session, id = "prgsBar_summaryOTUtablePer", value = (i/sampleNum)*100)
      } else {
        system(str_c(perTableOTUsComd1, samplename, perTableOTUsComd2, perTableOTUsComd3, 
                     samplename, perTableOTUsComd4, otuoutputFileName, perTableOTUsComd5))
        extractOTUsTable(str_c(parentnochimerasPathper, 'RST6_OTUsTablePerSample/'), otuoutputFileName)
        updateProgressBar(session, id = "prgsBar_summaryOTUtablePer", value = (i/sampleNum)*100)
      }
    }
  })
  
  
  # ==================================================================
  # Construct OTUs Table -- summary OTUs-matrix among all samples
  shinyDirChoose(input,
                 'rawnochimerasReadsAllPath',
                 roots = volumes)
  
  observeEvent(input$summaryOTUtableAll, {
    parentnochimerasPathall <- getParentFullPath(input$rawnochimerasReadsAllPath)
    rawnochimerassPathall <- getFullPath(input$rawnochimerasReadsAllPath)
    setwd(parentnochimerasPathall)
    sampleNames <- list.files(rawnochimerassPathall)
    sampleNum <- length(sampleNames)
    
    # creat & refresh RST7_OTUsTableAllSample
    if (dir.exists(str_c(parentnochimerasPathall,'RST7_OTUsTableAllSample/'))) {
      unlink(str_c(parentnochimerasPathall,'RST7_OTUsTableAllSample/'), recursive = TRUE)
    } else {
      dir.create(str_c(parentnochimerasPathall,'RST7_OTUsTableAllSample/'))
    }
    # creat & refresh RST7_AllNochimerasDereplicatedOTUs
    if (dir.exists(str_c(parentnochimerasPathall,'RST7_AllNochimerasDereplicatedOTUs/'))) {
      unlink(str_c(parentnochimerasPathall,'RST7_AllNochimerasDereplicatedOTUs/'), recursive = TRUE)
    } else {
      dir.create(str_c(parentnochimerasPathall,'RST7_AllNochimerasDereplicatedOTUs/'))
    }
    
    # extract all nochimerasFasta to one fastaFile
    system(str_c('python ./Scripts/ExtractAllNochimerasToOneFasta.py ', rawnochimerassPathall, 
                 ' ', parentnochimerasPathall,'RST7_AllNochimerasDereplicatedOTUs/'))
    updateProgressBar(session, id = "prgsBar_summaryNochimerasFastaAll", value = 50)
    
    # start to dereplicate AllNochimeras files
    dereplicateAllComd1 <- str_c('vsearch.exe --derep_fulllength ', 
                                 parentnochimerasPathall,'RST7_AllNochimerasDereplicatedOTUs/RST6_AllNochimeras.fasta')
    dereplicateAllComd2 <- str_c(' --sizeout --minuniquesize ', input$selectMinuniquesize)
    dereplicateAllComd3 <- str_c(' --output ', parentnochimerasPathall, 
                                 'RST7_AllNochimerasDereplicatedOTUs/RST6_AllDereplicated.fasta')
    system(str_c(dereplicateAllComd1, dereplicateAllComd2, dereplicateAllComd3))
    updateProgressBar(session, id = "prgsBar_summaryNochimerasFastaAll", value = 100)
    updateProgressBar(session, id = "prgsBar_summaryTotalOTUsFastaAll", value = 50)
    
    # start to cluster AllDereplicated files
    clusterAllComd1 <- str_c('vsearch.exe --cluster_fast ', 
                             parentnochimerasPathall,'RST7_AllNochimerasDereplicatedOTUs/RST6_AllDereplicated.fasta')
    clusterAllComd2 <- str_c(' --id ', input$selectClusterID)
    clusterAllComd3 <- str_c(' --relabel OTU_  --centroids ', 
                             parentnochimerasPathall,'RST7_AllNochimerasDereplicatedOTUs/RST6_AllOTUs.fasta')
    clusterAllComd4 <- str_c(' --iddef ', input$selectOTUsTableAllIDdef)
    if (input$selectOTUsTableAllIDdef == "ND") {
      system(str_c(clusterAllComd1, clusterAllComd2, clusterAllComd3))
      updateProgressBar(session, id = "prgsBar_summaryTotalOTUsFastaAll", value = 100)
    } else {
      system(str_c(clusterAllComd1, clusterAllComd2, clusterAllComd3, clusterAllComd4))
      updateProgressBar(session, id = "prgsBar_summaryTotalOTUsFastaAll", value = 100)
    }
    
    # construc OTUs table for each sample, reference -- RST6_AllOTUs.fasta
    allTableOTUsComd1 <- str_c('vsearch.exe --usearch_global ', rawnochimerassPathall)
    allTableOTUsComd2 <- str_c(' --id ', input$selectConstructOTUidAll)
    allTableOTUsComd3 <- str_c(' --db ', parentnochimerasPathall,'RST7_AllNochimerasDereplicatedOTUs/RST6_AllOTUs.fasta ')
    allTableOTUsComd4 <- str_c(' --otutabout ', parentnochimerasPathall, 'RST7_OTUsTableAllSample/')
    allTableOTUsComd5 <- str_c(' --iddef ', input$selectOTUsTableAllIDdef)
    for (i in 1:sampleNum) {
      samplename <- sampleNames[i]
      otuoutputFileName <- str_c(strsplit(samplename, '.fasta')[[1]], '.txt')
      if (input$selectOTUsTableAllIDdef == "ND") {
        system(str_c(allTableOTUsComd1, samplename, allTableOTUsComd2, 
                     allTableOTUsComd3, allTableOTUsComd4, otuoutputFileName))
        extractOTUsTable(str_c(parentnochimerasPathall, 'RST7_OTUsTableAllSample/'), otuoutputFileName)
        updateProgressBar(session, id = "prgsBar_summaryAllOTUtableEach", value = (i/sampleNum)*100)
      } else {
        system(str_c(allTableOTUsComd1, samplename, allTableOTUsComd2, 
                     allTableOTUsComd3, allTableOTUsComd4, otuoutputFileName, allTableOTUsComd5))
        extractOTUsTable(str_c(parentnochimerasPathall, 'RST7_OTUsTableAllSample/'), otuoutputFileName)
        updateProgressBar(session, id = "prgsBar_summaryAllOTUtableEach", value = (i/sampleNum)*100)
      }
    }
    
    # summary OTUsMatrix from all 'RST7_OTUsTableAllSample' files
    allOTUsNum <- length(getSequences(str_c(parentnochimerasPathall,'RST7_AllNochimerasDereplicatedOTUs/RST6_AllOTUs.fasta ')))
    OTUsEach.df <- data.frame(OTUs = str_c("OTU_", 1:allOTUsNum))
    for (i in 1:sampleNum) {
      samplename <- sampleNames[i]
      otuoutputFileName <- str_c(strsplit(samplename, '.fasta')[[1]], '.txt')
      tmpOTUsEach.df <- read.table( str_c(parentnochimerasPathall, 'RST7_OTUsTableAllSample/', otuoutputFileName) , 
                             header = F,sep = '\t',row.names = 1, stringsAsFactors = F)
      OTUs.DF <- data.frame(OTUs = rownames(tmpOTUsEach.df), Counter = rowSums(tmpOTUsEach.df))
      colnames(OTUs.DF) <- c("OTUs", strsplit(samplename, '.fasta')[[1]])
      OTUsEach.df <- merge(OTUsEach.df, OTUs.DF, by = "OTUs", all.x = TRUE)
      updateProgressBar(session, id = "prgsBar_summaryAllOTUsMatrix", value = (i/sampleNum)*100)
    }
    write.table(OTUsEach.df, str_c(parentnochimerasPathall, 
                                   'RST7_AllNochimerasDereplicatedOTUs/RST7_OTUsSamplesMatrix.txt'),
                quote = F, sep = '\t',row.names = F,col.names = T,na = "0")
    output$showAllOTUsMatrix <- DT::renderDataTable({ OTUsEach.df })
  })
  
  
  
  
  
  
  
  
  
  # ==================================================================
  observeEvent(input$assignTax2OTUperGlobal, {
    for (i in 1:10) {
      updateProgressBar(session, id = "prgsBar_assignTax2OTUperGlobal", value = i*10)
    }
  })
  
  
  
  # ==================================================================
  observeEvent(input$assignTax2OTUallGlobal, {
    for (i in 1:10) {
      updateProgressBar(session, id = "prgsBar_assignTax2OTUallGlobal", value = i*10)
    }
  })
  
  
  
  # ==================================================================
  observeEvent(input$assignTax2OTUperSintax, {
    for (i in 1:10) {
      updateProgressBar(session, id = "prgsBar_assignTax2OTUperSintax", value = i*10)
    }
  })
  
  
  
  
  # ==================================================================
  observeEvent(input$assignTax2OTUallSintax, {
    for (i in 1:10) {
      updateProgressBar(session, id = "prgsBar_assignTax2OTUallSintax", value = i*10)
    }
  })
  
  
  
  
  
  
  
  
}



shinyApp(ui = ui, server = server)