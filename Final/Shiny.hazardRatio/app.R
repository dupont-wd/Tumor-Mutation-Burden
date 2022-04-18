# shiny.hr.rcs.app.18.R

# Copyright (C) Vanderbilt University Medical Center

# See hr.cph.cancer.R for the non-shiny version of this program
# To publish this app it must be renamed app.R

# This program inputs values of CNAF and PMC from a Shiny app.
# It then calculates the hazard ratio and 95% confidence interval
# for a patient with these values of CNAF and PMC relative to
# those of a patient with CNAF = 0 and PMC = 1. It also draws a
# graph of the absoulte survival probability for a patient with
# with these values together with 95% confidence bands.

# The program also allows the cancer type and age to be entered.
# In this case no hazard ratios are calculated but the 15 year
# survival graph is drawn.

library(shiny)
library(survival)
library(rms)

pmc_cnaf_path <- readRDS("pmc_cnaf_path.RDS")
head(pmc_cnaf_path)
fate15  <- pmc_cnaf_path$fate15
table(fate15)
follow15 <- pmc_cnaf_path$follow15

# fate is an nx2 matrix giving the follow-up time and
# fate for each patient
fate <- Surv(follow15, fate15)
pmc_cnaf_path$fate <- Surv(follow15, fate15)
pmc <- pmc_cnaf_path$PointMutationCount
cnaf <- pmc_cnaf_path$CNAFractionGenomeAltered
cancer <- pmc_cnaf_path$cancer
pmc_cnaf_path$age <- pmc_cnaf_path$DiagnosisAge
stage <- pmc_cnaf_path$CancerStage
pmc_cnaf_path$stage <- pmc_cnaf_path$CancerStage
grade <- pmc_cnaf_path$TumorGrade
pmc_cnaf_path$grade <- pmc_cnaf_path$TumorGrade
head(pmc_cnaf_path$age)
logPM <- log(pmc)

# Embed logPM and cnaf in pmc_cnaf_path
pmc_cnaf_path$logPM <- logPM
pmc_cnaf_path$cnaf <- cnaf
 
# Estimate the hazard ratio associated with pmc_in and cnaf_in
# The ui defines the user interface.

#!!!!!!!!!!!!!!!!!!!! Start of the definition of ui !!!!!!!!!!!!!!!!!!!

ui <- fluidPage(

    titlePanel(h1("Hazard ratio and survival calculator based on SM count, CNA fraction, age, stage and grade"
        ,style ="font-size:150%; font-weight: bold", align = "center")),

# The sidebarPanel consits of all input and output on the left of
# the shiny app.

    sidebarPanel(
    
# Copy the line below to make a number input box into the UI.
# input$pmc_in and input$cnaf_in are the values of PMC and CNAF
# entered by the user.

        numericInput("pmc_in",  label = h3("SM count", style ="font-size:100%; font-weight: bold"), value = 1, min = 0),
        numericInput("cnaf_in",  label = h3("CNA fraction", style ="font-size:100%; font-weight: bold"), 
            value = 0, min = 0, max = 1, step =.01 ),
        selectInput("cancer_in", "Cancer", 
            list("Not chosen",
                "Adrenocortical Carcinoma",
                "Bladder Cancer",
                "Breast Cancer",
                "Cervical Squamous Cell Carcinoma",
                "Cholangiocarcinoma",
                "Colorectal Adenocarcinoma",
                "Diffuse Large B-Cell Lymphoma",
                "Esophageal Adenocarcinoma",
                "Glioblastoma Multiforme",
                "Head and Neck Squamous Cell Carcinoma",
                "Kidney Chromophobe",
                "Kidney Renal Clear Cell Carcinoma",
                "Kidney Renal Papillary Cell Carcinoma",
                "Acute Myeloid Leukemia",
                "Brain Lower Grade Glioma",
                "Liver Hepatocellular Carcinoma",
                "Lung Adenocarcinoma",
                "Lung Squamous Cell Carcinoma",
                "Mesothelioma",
                "Ovarian Serous Cystadenocarcinoma",
                "Pancreatic Adenocarcinoma",
                "Pheochromocytoma and Paraganglioma",
                "Prostate Adenocarcinoma",
                "Sarcoma",
                "Skin Cutaneous Melanoma",
                "Stomach Adenocarcinoma",
                "Testicular Germ Cell Tumors",
                "Thyroid Carcinoma",
                "Thymoma",
                "Uterine Corpus Endometrial Carcinoma",
                "Uterine Carcinosarcoma",
                "Uveal Melanoma")
        ),
         numericInput("age_in", label = h3("Age", style ="font-size:100%; font-weight: bold" ), value = " ", min = 1),
        selectInput("stage_in", "Stage", 
           list("Not chosen",
                "1",
                "2",
                "3",
                "4")
        ),
        selectInput("grade_in", "Grade", 
           list("Not chosen",
                "1",
                "2",
                "3",
                "4")
        ),

        fluidRow(column(3, verbatimTextOutput("value"))),
            
# The output order is as follows. 
# The order of the textOutput commands will determine the order that
# this output is given.

        textOutput("badInput1"),
        textOutput("badInput2"),
        textOutput("badInput3"),
        textOutput("badInput4"),
        textOutput("badInput5"),      
        textOutput("hrtext"),
            
# hr, calculated by output$hr below will be the derived hazard ratio
        textOutput("hr"), 
        textOutput("lbtext"),
            
# lb will be the lower bound of the 95% CI calculated by output$lb   
        textOutput("lb"),
        textOutput("ubtext"),
            
# ub will be the upper bound of this inteval   
        textOutput("ub"),
        textOutput("result"),
        textOutput("numresult"),
            
# inprogress is a temporary field
        textOutput("inprogress"),
        
# line draws a line after the output
        textOutput("line"),
        
# url gives a url reference
        uiOutput("url")
    ),
  

# The mainPanel is where Shiny draws graphs  
    mainPanel(
        plotOutput('plot')
    )
) 

# server is the server function where calculations are performed.
# ui and server are combined together in the shinyApp command in
# the last line of this program.

#!!!!!!!!!!!!!!!!!! Start of the definition of server !!!!!!!!!!!!!!!!!!!!!!!!
server <- function(input, output, session) {

    
# badInput is an error message that is displayed when the inputs are invalid

    output$badInput1 <- renderText({
        if (!(isTruthy(input$pmc_in)) | !(isTruthy(input$cnaf_in))) {
            badInput1 <- "PMC must be >= 1 & CNAF must be > 0 & <= 1"
        }    
        else {    
            if (!(input$pmc_in >= 1  & input$cnaf_in >= 0 & input$cnaf_in <= 1)) {
                badInput1 <- "PMC must be >= 1 & CNAF must be > 0 & <= 1"
            }
        }
    })
    output$badInput2 <- renderText({
        if ((isTruthy(input$age_in)) & input$cancer_in == "Not chosen") {
            badInput2 <-  "Age must be missing when cancer is not chosen "  
        }
    })
    output$badInput4 <- renderText({
        if (input$stage_in != "Not chosen" & input$cancer_in == "Not chosen") {
            badInput2 <-  "Stage can only be chosen when cancer is chosen "  
        }
    })
    output$badInput5 <- renderText({
        if ((input$grade_in != "Not chosen" & input$cancer_in == "Not chosen") 
          | (input$grade_in != "Not chosen" & input$stage_in == "Not chosen")) {
            badInput2 <-  "grade can only be chosen when cancer and stage are chosen "  
        }
    })
    
# We only want to calculate hazard ratios for valid inputs when no cancer is specified
# and age is missing. 

    output$hrtext <- renderText({
        if ((isTruthy(input$pmc_in)) & (isTruthy(input$cnaf_in))) {
            if (input$pmc_in >= 1 & input$cnaf_in >= 0 & input$cnaf_in <= 1 & 
                input$cancer_in == "Not chosen" & input$stage_in == "Not chosen"  
                & input$grade_in == "Not chosen" & !isTruthy(input$age_in) ) {
                hrtext <- "Mortal hazard relative to patients with SM count = 1 and CNA fraction = 0 "
            } 
        }
    })

    output$hr <- renderText({
        if ((isTruthy(input$pmc_in)) & (isTruthy(input$cnaf_in))) { 
            if (input$pmc_in < 1 | input$cnaf_in < 0 | input$cnaf_in >= 1  ) { 
                hr <- " "
            }   
            else {
                if (input$cancer_in == "Not chosen" & input$stage_in == "Not chosen" 
                    & input$grade_in == "Not chosen" & !isTruthy(input$age_in)) {
            
############## Start of code to calculate hazard ratios when no cancer or stage are chosen ##########            
    
# Regress fate against a restricted cubic 
# spline model with 3 knots for log(PMC) and 4 knots for CNAF. The 
# knot locations are given explicity and are the default Harrel knots as 
# calcualted by Stata. These are identical to the cph knot locations to
# 7 significant figures
 
                    model <- cph(formula =fate~rcs(logPM,c( 2.564949,   4.077538,   5.897154 )) 
                        + rcs(cnaf, c( .0001,     .11968,     .29422,     .68018 )),data = pmc_cnaf_path)

# Calculate the hazard for someone with pmc =input$pmc_in & cnaf = input$cnaf_in relative 
# to someone with pmc =1 and cnaf= 0. model refers to the hazard regression given above.
# loghr is a list giving the log hazard ratios, confidence intervals and other statistics for
# each patient. The contrast command has two lists of inputs. The first and second list
# give the numerator and denominator of the hazard ratio, respectively.

                    cnaf_in <- input$cnaf_in
                    if(cnaf_in ==0) {
                        cnaf_in <- 0.0001
                    }    
                    pmc_in <- input$pmc_in
                    if(pmc_in ==1) {
                        pmc_in <- 1.0001
                    }    
                    loghr <- contrast(model,list(logPM = log(pmc_in), cnaf = cnaf_in),list(logPM = 0, cnaf=0))
 
# loghr is a list. To extract an element from a list you need to use double brackets
# as shown below. The log hazard ratio is given in by the first element of this list. The 
# lower and upper bound of the 95% CI are given in positions 3 and 4, respectively.

                    hrval <- signif(exp(loghr[[1]]),digits=3)
                    hr <- hrval
                    hr
                }
            } 
        }
    })
    output$lbtext <- renderText({
        if ((isTruthy(input$pmc_in)) & (isTruthy(input$cnaf_in))) { 
            if (input$cancer_in == "Not chosen"  & input$stage_in == "Not chosen" 
                & input$grade_in == "Not chosen" & !isTruthy(input$age_in)) {
                if ((input$pmc_in >= 1 & input$cnaf_in >= 0 & input$cnaf_in <= 1)) {
                    lbtext <- "Lower bound of 95% CI"
                } 
            } 
        }
    })
    
    output$lb <- renderText({
        if ((isTruthy(input$pmc_in)) & (isTruthy(input$cnaf_in))) {     
            if (input$cancer_in == "Not chosen" & input$stage_in == "Not chosen"  
                & input$grade_in == "Not chosen"& !isTruthy(input$age_in)) {
                if ((input$pmc_in < 1 | input$cnaf_in < 0 | input$cnaf_in >= 1)) {
                    lb <- " "
                } 
                else { 
                    cnaf_in <- input$cnaf_in
                    if(cnaf_in ==0) {
                        cnaf_in <- 0.0001
                    } 
                    pmc_in <- input$pmc_in
        	        if(pmc_in ==1) {
        	            pmc_in <- 1.0001
        	        }    
                    model <- cph(formula =fate~rcs(logPM,c(2.564949,   4.077538,   5.897154)) 
                        + rcs(cnaf, c(.0001,     .11968,     .29422,     .68018)),data = pmc_cnaf_path)
                    loghr <- contrast(model,list(logPM = log(pmc_in), cnaf=cnaf_in),list(logPM = 0, cnaf=0))
                    hrval <- signif(exp(loghr[[1]]),digits=3)
                    lbval <-  signif(exp(loghr[[3]]),digits=2)
                    if (lbval > hrval) {
                        lbval <- hrval
                    }
                    lb <- lbval
                    lb
                }
            }
        }
    })
    output$ubtext <- renderText({
        if ((isTruthy(input$pmc_in)) & (isTruthy(input$cnaf_in))) {     
            if (input$cancer_in == "Not chosen" & input$stage_in == "Not chosen" 
                & input$grade_in == "Not chosen" & !isTruthy(input$age_in)) { 
                if ((input$pmc_in >= 1 & input$cnaf_in >= 0 & input$cnaf_in <= 1)) {
                    ubtext <- "Upper bound of 95% CI"
                }
            }
        }    
    })
      
    output$ub <- renderText({ 
        if ((isTruthy(input$pmc_in)) & (isTruthy(input$cnaf_in))) {         
            if (input$cancer_in == "Not chosen" & input$stage_in == "Not chosen" 
                & input$grade_in == "Not chosen" & !isTruthy(input$age_in)) {
                if ((input$pmc_in < 1 | input$cnaf_in < 0 | input$cnaf_in >= 1)) {
                    ub <- " "
                } 
                else {   
                    cnaf_in <- input$cnaf_in
                    if(cnaf_in ==0) {
                        cnaf_in <- 0.0001
                    }
                    pmc_in <- input$pmc_in
                    if(pmc_in ==1) {
                       pmc_in <- 1.0001
                    }    
                    model <- cph(formula =fate~rcs(logPM,c(2.564949,   4.077538,   5.897154)) 
                        + rcs(cnaf, c(.0001,     .11968,     .29422,     .68018)),data = pmc_cnaf_path)
                    loghr <- contrast(model,list(logPM = log(pmc_in), cnaf=cnaf_in),list(logPM = 0, cnaf=0))
                    hrval <- signif(exp(loghr[[1]]),digits=3)
                    ubval <-  signif(exp(loghr[[4]]),digits=2)
                    if (ubval < hrval) {
                        ubval <- hrval
                    }
                    ub <- ubval
                    ub
                }
            } 
        }
    })     
    output$line <-renderText({
        line <- "_______________________________ "
    })
    url <- a("here", href="https://biostat.app.vumc.org/wiki/Main/TumorMutationBurden")
    output$url <- renderUI({
      tagList("Click ", url, " for more information.")
    })
         
    #!!!!!!!!!!!!!!!!!! Start of survival plots code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    output$plot <- renderPlot({
    # Check that TMB covariates are given and are valid
        if ((isTruthy(input$pmc_in)) & (isTruthy(input$cnaf_in)) &
            input$pmc_in >= 1 & input$cnaf_in >= 0 &  input$cnaf_in <= 1 ) {    
            
   #@@@@@@@@@@ Plot survival curve when TMB covariates are valid, age is missing and cancer, stage and grade are not chosen @@@@        
            if (input$cancer_in == "Not chosen" & !isTruthy(input$age_in)) { 
                if (input$stage_in == "Not chosen"  & input$grade_in == "Not chosen") {
                    cnaf_in <- input$cnaf_in
                    if(cnaf_in ==0) {
                        cnaf_in <- 0.0001
                    }
                    pmc_in <- input$pmc_in
                    if(pmc_in ==1) {
                        pmc_in <- 1.0001
                    }    
                    model <- cph(formula =fate~rcs(logPM,c(2.564949,   4.077538,   5.897154)) 
                        + rcs(cnaf, c(.0001,     .11968,     .29422,     .68018)) , data = pmc_cnaf_path, x=TRUE, y=TRUE)
                    p <- survplot(model, logPM = log(input$pmc_in), cnaf = input$cnaf_in, conf.int=.95,col='red', col.fill='yellow', 
                        xlab="Years since diagnosis", cex.xlab=1.3,  time.inc = 1, 
                        ylab="Probability of overall survival", lwd=3, grid=TRUE)
                    print(p) 
                }
            }
   #@@@@@@@@@ End of this plot
            
            else {
   # First, define labels for cancer, stage and grade        
                if (input$cancer_in != "Not chosen" ) {
                    cancer.labels = c(
                        "Adrenocortical Carcinoma",
                        "Bladder Cancer",
                        "Breast Cancer",
                        "Cervical Squamous Cell Carcinoma",
                        "Cholangiocarcinoma",
                        "Colorectal Adenocarcinoma",
                        "Diffuse Large B-Cell Lymphoma",
                        "Esophageal Adenocarcinoma",
                        "Glioblastoma Multiforme",
                        "Head and Neck Squamous Cell Carcinoma",
                        "Kidney Chromophobe",
                        "Kidney Renal Clear Cell Carcinoma",
                        "Kidney Renal Papillary Cell Carcinoma",
                        "Acute Myeloid Leukemia",
                        "Brain Lower Grade Glioma",
                        "Liver Hepatocellular Carcinoma",
                        "Lung Adenocarcinoma",
                        "Lung Squamous Cell Carcinoma",
                        "Mesothelioma",
                        "Ovarian Serous Cystadenocarcinoma",
                        "Pancreatic Adenocarcinoma",
                        "Pheochromocytoma and Paraganglioma",
                        "Prostate Adenocarcinoma",
                        "Sarcoma",
                        "Skin Cutaneous Melanoma",
                        "Stomach Adenocarcinoma",
                        "Testicular Germ Cell Tumors",
                        "Thyroid Carcinoma",
                        "Thymoma",
                        "Uterine Corpus Endometrial Carcinoma",
                        "Uterine Carcinosarcoma",
                        "Uveal Melanoma")
                    cancer.f <- factor(cancer, labels = cancer.labels )
                    cancer_in <-input$cancer_in
                    pmc_cnaf_path$cancer.f <- factor(pmc_cnaf_path$cancer, labels = cancer.labels )                 
                }
                if (input$stage_in != "Not chosen") {
                    stage.labels = c(
                        "1",
                        "2",
                        "3",
                        "4")
                    stage.f <- factor(stage, labels = stage.labels )
                    stage_in <-input$stage_in
                    pmc_cnaf_path$stage.f <- factor(pmc_cnaf_path$stage, labels = stage.labels )                                      
                }
                if (input$grade_in != "Not chosen") {
                    grade.labels = c(
                        "1",
                        "2",
                        "3",
                        "4")
                    grade.f <- factor(grade, labels = grade.labels )
                    grade_in <-input$grade_in
                    pmc_cnaf_path$grade.f <- factor(pmc_cnaf_path$grade, labels = grade.labels )                                      
                }
                
#@@@@@@@@@@ Plot survival curve when TMB covariates are valid, cancer is chosen, age is missing and stage and grade are not chosen @@@@                  
                if (input$cancer_in != "Not chosen" & !isTruthy(input$age_in) 
                    & input$stage_in == "Not chosen" & input$grade_in == "Not chosen") {
                    cnaf_in <- input$cnaf_in
                    if(cnaf_in ==0) {
                        cnaf_in <- 0.0001
                    }
                    pmc_in <- input$pmc_in
                    if(pmc_in ==1) {
                        pmc_in <- 1.0001
                    }  
    
# Regress fate against logPM cnaf and cancer                
                    model <- cph(formula =fate~rcs(logPM,c(2.564949,   4.077538,   5.897154)) 
                        + rcs(cnaf, c(.0001,     .11968,     .29422,     .68018))
                        + cancer.f, data = pmc_cnaf_path, x=TRUE, y=TRUE)
                    p <- survplot(model, logPM = log(pmc_in), cnaf = cnaf_in, cancer.f=cancer_in,
                             conf.int=.95,col='red', col.fill='yellow', 
                             xlab="Years since diagnosis",  cex.xlab=1.3, time.inc = 1, 
                             ylab="Probability of overall survival", lwd=3, grid=TRUE)
                    print(p) 
                }
#@@@@@@@@@@ Plot survival curve when TMB covariates are valid, cancer and stage are chosen, age is missing and grade is not chosen @@@@                  
                
                if (input$cancer_in != "Not chosen" & input$stage_in != "Not chosen" 
                    & input$grade_in == "Not chosen" &  !isTruthy(input$age_in)) {                     

# tapply creates a list of all the cancer labels in pmc_cnaf_path$cancer.f. For each cancer label in this list it creates
# a sublist of all the distinct stages in pmc_cnaf_path$stage that are found in patients with this cancer. 
# It then applies some function to this list of sublists. In this example the function, called FUN, is
# {all(is.na(x))}. It returns a vector with one value for every element of pmc_cnaf_path$cancer.f. Each
# value in this vector is TRUE if all values of pmc_cnaf_path$stage are missing for the specified cancer and is false otherwise.
# select_labels negates this execution of tapply

                    select_labels <-!tapply(pmc_cnaf_path$stage, pmc_cnaf_path$cancer.f, FUN=function(x) {all(is.na(x))})
                    select_labels
                    
# Redefine cancer.include to be those labels in cancer.labels for which at least one record associated with
# the cancer has a non-missing stage.

                    cancer.include <- cancer.labels[select_labels]
  
# Restrict the indicator covariates in cancer.f to those covariates in cancer.include
                    data <- subset(pmc_cnaf_path, cancer.f %in% cancer.include, select=c(logPM,cnaf,fate,stage.f,cancer.f)) 
                    data$cancer.f <- droplevels(data$cancer.f)

# calculate survival curve for someone with a specified cancer-type, TMB and stage
# Regress fate against logPM cancer.f and cnaf
                    print(model <- cph(formula =fate~rcs(logPM,c(2.6390574, 4.0775375, 5.8998976)) 
                        + rcs(cnaf, c(0, .1165, .2924, .678)) + stage.f+ cancer.f,
                    data = data,surv=TRUE, x=TRUE, y=TRUE))
                    mycancer =  input$cancer_in
# mycancer will be true if mycancer is in the cancer list cancer.include                    
                    mycancer %in% cancer.include
                    cnaf_in <- input$cnaf_in
		        if(cnaf_in ==0) {
		            cnaf_in <- 0.0001
		        }
		        pmc_in <- input$pmc_in
		            if(pmc_in ==1) {
		                pmc_in <- 1.0001
		            }
# We only want to draw a survival plot if there are non-missing values of stage for the cancer of interest		            
                     if (mycancer %in% cancer.include) {
                        p2 <- survplot(model, logPM = log(pmc_in), cnaf = cnaf_in, stage.f=input$stage_in, cancer.f=mycancer 
                            , conf.int=.95,col='red', col.fill='yellow', 
                            xlab="Years since diagnosis",  cex.xlab=1.3, time.inc = 1, 
                            ylab="Probability of overall survival", lwd=3, grid=TRUE)
                        print(p2)
                    output$badInput3 <- renderText({
                        badInput3 <-  ""
                    })
                    }                       
                    else {
                        output$badInput3 <- renderText({
                            if (!(mycancer %in% cancer.include) &  isTruthy(input$stage_in)) {
                                badInput3 <-  "Pan-cancer Atlas does not contain the data needed for this survival curve"
                            }    
                        })
                    }    
                }

#@@@@@@@@@@ Plot survival curve when TMB covariates are valid, cancer and age are specified and stage and grade are not @@@@                  

                if (input$cancer_in != "Not chosen" &  isTruthy(input$age_in) 
                    & input$stage_in == "Not chosen" & input$grade_in == "Not chosen" ) {
                        
# Age is completely missing for Acute Myeloid Leukemia and Uveal Melanoma. We must
# eliminate records with these cancers from the analysis and eliminate these
# cancers from the list of cancers in the model

# tapply creates a list of all the cancer labels in pmc_cnaf_path$cancer.f. For each cancer label in this list it creates
# a sublist of all the distinct ages in pmc_cnaf_path$age that are found in patients with this cancer. 
# It then applies some function to this list of sublists. In this example the function, called FUN, is
# {all(is.na(x))}. It returns a vector with one value for every element of pmc_cnaf_path$cancer.f. Each
# value in this vector is TRUE if all values of pmc_cnaf_path$age are missing for the specified cancer and is false otherwise.
# select_labels negates this execution of tapply
                    select_labels <-!tapply(pmc_cnaf_path$age, pmc_cnaf_path$cancer.f, FUN=function(x) {all(is.na(x))})
                    select_labels
# Redefine cancer.include to be those labels in cancer.labels for which at least one record associated with
# the cancer has a non-missing age.
                    cancer.include <- cancer.labels[select_labels]
  
# Restrict the indicator covariates in cancer.f to those covariates in cancer.include
                    data <- subset(pmc_cnaf_path, cancer.f %in% cancer.include, select=c(logPM,cnaf,fate,age,cancer.f))
                    data$cancer.f <- droplevels(data$cancer.f)

# calculate survival curve for someone with a specified cancer-type, TMB and age
# Regress fate against logPM cancer.f and cnaf
                    print(model <- cph(formula =fate~rcs(logPM,c(2.6390574, 4.0775375, 5.8998976)) 
                        + rcs(cnaf, c(0, .1165, .2924, .678)) + age+ cancer.f,
                    data = data,surv=TRUE, x=TRUE, y=TRUE))
                    mycancer =  input$cancer_in
# mycancer will be true if mycancer is in the cancer list cancer.include                    
                    mycancer %in% cancer.include
                    cnaf_in <- input$cnaf_in
		        if(cnaf_in ==0) {
		            cnaf_in <- 0.0001
		        }
		        pmc_in <- input$pmc_in
		            if(pmc_in ==1) {
		                pmc_in <- 1.0001
		            }
# We only want to draw a survival plot if there are non-missing values of age for the cancer of interest		            
                     if (mycancer %in% cancer.include) {
                        p2 <- survplot(model, logPM = log(pmc_in), cnaf = cnaf_in, cancer.f=mycancer, age=input$age_in 
                            , conf.int=.95,col='red', col.fill='yellow', 
                            xlab="Years since diagnosis",  cex.xlab=1.3, time.inc = 1, 
                            ylab="Probability of overall survival", lwd=3, grid=TRUE)
                        print(p2)
                    output$badInput3 <- renderText({
                        badInput3 <-  ""
                    })
                    }                       
                    else {
                        output$badInput3 <- renderText({
                            if (!(mycancer %in% cancer.include) &  isTruthy(input$age_in)) {
                                badInput3 <-  "Pan-cancer Atlas does not contain the data needed for this survival curve"
                            }    
                        })
                    }    
                }
#@@@@@@@@@@ Plot survival curve when TMB covariates are valid, cancer, stage and age are given and grade is not @@@@  

                if (input$cancer_in != "Not chosen" & input$stage_in != "Not chosen" 
                    &  isTruthy(input$age_in) & input$grade_in == "Not chosen") {
                        
# tapply creates a list of all the cancer labels in pmc_cnaf_path$cancer.f. For each cancer label in this list it creates
# a sublist of all the distinct stages in pmc_cnaf_path$stage that are found in patients with this cancer. 
# It then applies some function to this list of sublists. In this example the function, called FUN, is
# {all(is.na(x))}. It returns a vector with one value for every element of pmc_cnaf_path$cancer.f. Each
# value in this vector is TRUE if all values of pmc_cnaf_path$stage are missing for the specified cancer and is false otherwise.
# select_labels negates this execution of tapply

                    select_labels <-!tapply(pmc_cnaf_path$stage, pmc_cnaf_path$cancer.f, FUN=function(x) {all(is.na(x))})
                    select_labels
                    
# Redefine cancer.include to be those labels in cancer.labels for which at least one record associated with
# the cancer has a non-missing stage.

                    cancer.include <- cancer.labels[select_labels]
  
# Restrict the indicator covariates in cancer.f to those covariates in cancer.include
                    data <- subset(pmc_cnaf_path, cancer.f %in% cancer.include, select=c(logPM,cnaf,fate,age,stage.f,cancer.f)) 
                    data$cancer.f <- droplevels(data$cancer.f)

# Regress fate against logPM, cnaf, age stage.f and cancer.f 
                    print(model <- cph(formula =fate~rcs(logPM,c(2.6390574, 4.0775375, 5.8998976)) 
                        + rcs(cnaf, c(0, .1165, .2924, .678)) + age + stage.f+ cancer.f,
                    data = data,surv=TRUE, x=TRUE, y=TRUE))
                    mycancer =  input$cancer_in
                    mystage = input$stage_in
                    myage = input$age_in
# mycancer will be true if mycancer is in the cancer list cancer.include                    
                    mycancer %in% cancer.include
                    cnaf_in <- input$cnaf_in
		        if(cnaf_in ==0) {
		            cnaf_in <- 0.0001
		        }
		        pmc_in <- input$pmc_in
		            if(pmc_in ==1) {
		                pmc_in <- 1.0001
		            }
# We only want to draw a survival plot if there are non-missing values of age for the cancer of interest		            
                     if (mycancer %in% cancer.include) {
                        p2 <- survplot(model, logPM = log(pmc_in), cnaf = cnaf_in, age=myage,  
                            stage.f=mystage, cancer.f=mycancer, conf.int=.95,col='red', col.fill='yellow', 
                            xlab="Years since diagnosis",  cex.xlab=1.3, time.inc = 1, 
                            ylab="Probability of overall survival", lwd=3, grid=TRUE)
                        print(p2)
                    output$badInput3 <- renderText({
                        badInput3 <-  ""
                    })
                    }                       
                    else {
                        output$badInput3 <- renderText({
                            if (!(mycancer %in% cancer.include) ) {
                                badInput3 <-  "Pan-cancer Atlas does not contain the data needed for this survival curve"
                            }    
                        })
                    }    
                }
#@@@@@@@@@@ Plot survival curve when TMB covariates are valid, cancer, stage and grade are given but age is not @@@@  

                if (input$cancer_in != "Not chosen" & input$stage_in != "Not chosen" 
                    &  !isTruthy(input$age_in) & input$grade_in != "Not chosen") {
                        
# Create a variable pmc_cnaf_path$stage_grade that is missing whenever stage or  grade is missing
# stage_grade <- subset(pmc_cnaf_path, !is.na(stage) & !is.na(grade) , select= c( cancer.f))

# split creates a list of all the cancer labels in pmc_cnaf_path$cancer.f. For each cancer label in this list it creates
# a sublist of all the records that are found in patients with this cancer.
# sapply is passed the list that split created.
# It then applies some function to this list of sublists. In this example the function, passed to the FUN argument, is
# {all(is.na(x))}. It returns a vector with one value for every element of pmc_cnaf_path$cancer.f. Each
# value in this vector is TRUE if all values of either pmc_cnaf_path$stage or pmc_cnaf_path$grade are missing for the 
# specified cancer and is false otherwise.
# select_labels stores the complement of this result produced by sapply.

                    select_labels <- !sapply(split(pmc_cnaf_path[,c("stage","grade")], pmc_cnaf_path$cancer.f), 
                                             FUN=function(x) { all(is.na(x$stage) | is.na(x$grade))})
                    print(select_labels)
                    
# Redefine cancer.include to be those labels in cancer.labels for which at least one record associated with
# the cancer has a non-missing stage.

                    cancer.include <- cancer.labels[select_labels]
                    cancer.include
  
# Restrict the indicator covariates in cancer.f to those covariates in cancer.include
                    data <- subset(pmc_cnaf_path, cancer.f %in% cancer.include, 
                             select=c(logPM,cnaf,fate,age,stage.f,grade.f,cancer.f)) 
                    data$cancer.f <- droplevels(data$cancer.f)
                    data$stage.f <- droplevels(data$stage.f)
                    data$grade.f <- droplevels(data$grade.f)

# Regress fate against logPM, cnaf, age stage.f grade.f and cancer.f 
                    print(model <- cph(formula =fate~rcs(logPM,c(2.6390574, 4.0775375, 5.8998976)) 
                        + rcs(cnaf, c(0, .1165, .2924, .678))  + stage.f+ grade.f + cancer.f,
                    data = data,surv=TRUE, x=TRUE, y=TRUE))
                    mycancer =  input$cancer_in
                    mystage = input$stage_in
                    myage = input$age_in
                    mygrade = input$grade_in
# mycancer will be true if mycancer is in the cancer list cancer.include                    
                    mycancer %in% cancer.include
                    cnaf_in <- input$cnaf_in
		        if(cnaf_in ==0) {
		            cnaf_in <- 0.0001
		        }
		        pmc_in <- input$pmc_in
		            if(pmc_in ==1) {
		                pmc_in <- 1.0001
		            }
# We only want to draw a survival plot if there are non-missing values of stage and grade for the cancer of interest		            
                     if (mycancer %in% cancer.include) {
                        p2 <- survplot(model, logPM = log(pmc_in), cnaf = cnaf_in,   
                            stage.f=mystage, cancer.f=mycancer, grade.f = mygrade,
                            conf.int=.95,col='red', col.fill='yellow', 
                            xlab="Years since diagnosis",  cex.xlab=1.3, time.inc = 1, 
                            ylab="Probability of overall survival", lwd=3, grid=TRUE)
                        print(p2)
                    output$badInput3 <- renderText({
                        badInput3 <-  ""
                    })
                    }                       
                    else {
                        output$badInput3 <- renderText({
                            if (!(mycancer %in% cancer.include) ) {
                                badInput3 <-  "Pan-cancer Atlas does not contain the data needed for this survival curve"
                            }    
                        })
                    }    
                }
#@@@@@@@@@@ Plot survival curve when TMB covariates are valid, cancer, stage, grade and age are given. @@@@  
                if (input$cancer_in != "Not chosen" & input$stage_in != "Not chosen" 
                    &  isTruthy(input$age_in) & input$grade_in != "Not chosen") {
                        
# Create a variable pmc_cnaf_path$stage_grade that is missing whenever stage, grade or age3 is missing
# stage_grade <- subset(pmc_cnaf_path, !is.na(stage) & !is.na(grade) , select= c( cancer.f))

# split creates a list of all the cancer labels in pmc_cnaf_path$cancer.f. For each cancer label in this list it creates
# a sublist of all the records that are found in patients with this cancer.
# sapply is passed the list that split created.
# It then applies some function to this list of sublists. In this example the function, passed to the FUN argument, is
# {all(is.na(x))}. It returns a vector with one value for every element of pmc_cnaf_path$cancer.f. Each
# value in this vector is TRUE if all values of either pmc_cnaf_path$stage, pmc_cnaf_path$grade or pmc_cnaf_path$age 
# are missing for the specified cancer and is false otherwise.
# select_labels stores the complement of this result produced by sapply.

                    select_labels <- !sapply(split(pmc_cnaf_path[,c("stage","grade", "age")], pmc_cnaf_path$cancer.f), 
                                             FUN=function(x) { all(is.na(x$stage) | is.na(x$grade))})
                    print(select_labels)
                    
# Redefine cancer.include to be those labels in cancer.labels for which at least one record associated with
# the cancer has a non-missing stage.

                    cancer.include <- cancer.labels[select_labels]
                    cancer.include
  
# Restrict the indicator covariates in cancer.f to those covariates in cancer.include
                    data <- subset(pmc_cnaf_path, cancer.f %in% cancer.include, 
                             select=c(logPM,cnaf,fate,age,stage.f,grade.f,cancer.f)) 
                    data$cancer.f <- droplevels(data$cancer.f)
                    data$stage.f <- droplevels(data$stage.f)
                    data$grade.f <- droplevels(data$grade.f)
                    

# Regress fate against logPM, cnaf, age stage.f grade.f and cancer.f 
                    print(model <- cph(formula =fate~rcs(logPM,c(2.6390574, 4.0775375, 5.8998976)) 
                        + rcs(cnaf, c(0, .1165, .2924, .678)) +age  + stage.f+ grade.f + cancer.f,
                    data = data,surv=TRUE, x=TRUE, y=TRUE))
                    mycancer =  input$cancer_in
                    mystage = input$stage_in
                    myage = input$age_in
                    mygrade = input$grade_in
# mycancer will be true if mycancer is in the cancer list cancer.include                    
                    mycancer %in% cancer.include
                    cnaf_in <- input$cnaf_in
		        if(cnaf_in ==0) {
		            cnaf_in <- 0.0001
		        }
		        pmc_in <- input$pmc_in
		            if(pmc_in ==1) {
		                pmc_in <- 1.0001
		            }
# We only want to draw a survival plot if there are non-missing values of stage and grade for the cancer of interest		            
                     if (mycancer %in% cancer.include) {
                        p2 <- survplot(model, logPM = log(pmc_in), cnaf = cnaf_in,   
                            stage.f=mystage, cancer.f=mycancer, grade.f = mygrade, age=myage,
                            conf.int=.95,col='red', col.fill='yellow', 
                            xlab="Years since diagnosis",  cex.xlab=1.3, time.inc = 1, 
                            ylab="Probability of overall survival", lwd=3, grid=TRUE)
                        print(p2)
                    output$badInput3 <- renderText({
                        badInput3 <-  ""
                    })
                    }                       
                    else {
                        output$badInput3 <- renderText({
                            if (!(mycancer %in% cancer.include) ) {
                                badInput3 <-  "Pan-cancer Atlas does not contain the data needed for this survival curve"
                            }    
                        })
                    }    
                }
            }
        }
    } ,height=600)
}
    
shinyApp(ui, server)
