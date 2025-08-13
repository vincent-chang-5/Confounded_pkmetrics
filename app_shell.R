library(mrgsolve)
library(shiny)
library(ggplot2)
library(dplyr)
library(survival)
library(survminer)
library(tidyverse)
library(bslib)

model_list <- list.files(modlib(), pattern = "^pk.*\\.cpp$")
model_list <- gsub("\\.cpp$", "", model_list)

ui <- fluidPage(
  titlePanel("Time-Varying PK Metrics Bias"),
  sidebarLayout(
    sidebarPanel(
      numericInput("cl_choice", "Select a Clearance:", value = 20, min = 1, max = 100),
      numericInput("ii_choice", "Select Dosing Interval (days):", value = 21, min = 1, max = 100),
      sliderInput("sim_end", "Select a Simulation End Time (days):", value = 365, min = 1, max = 1000),
      radioButtons("pkmetric", "PK metric selection:",
                   c("CMAX" = "CMAX",
                     "CMIN" = "CMIN",
                     "CAVG" = "CAVG")),
      plotOutput("sim_plot"),
    ),
    mainPanel(navset_tab(
      nav_panel("Null Cox Model", plotOutput("km_plot_null"),
                plotOutput("km_plot_nullq"),
                plotOutput("km_plot_nullqt")),
      nav_panel("ER Cox Model", plotOutput("km_plot_corr"),
                plotOutput("km_plot_corrq"),
                plotOutput("km_plot_corrqt")),
      nav_panel("Null LR Model", plotOutput("lr_plot_nullq"),
                plotOutput("lr_plot_nullqt")),
      nav_panel("ER LR Model", plotOutput("lr_plot_corrq"),
                plotOutput("lr_plot_corrqt"))
    ))
  )
)

server <- function(input, output) {

  sim_data <- reactive({
    mod <- mread("model.cpp")

    # Adjust input$cl_choice for selected model
    updated_mod <- param(mod, TVCL = input$cl_choice)

    # Adjust input$sim_end for simulation end time
    sim_data <- updated_mod %>%
      ev(amt = 100, ii = 24*input$ii_choice, addl = 50) %>%
      mrgsim_df(nid = 500, end = 24*input$sim_end, delta = 1)
  })

  output$sim_plot <- renderPlot({
    ggplot(sim_data() %>%
             select(ID, time, IPRED, CAVG) %>%
             pivot_longer(c("IPRED","CAVG"), names_to = "pk", values_to = "value") %>%
             group_by(time,pk) %>%
             summarise(median = median(value),
                       upper = quantile(value, probs = 0.975),
                       lower = quantile(value, probs = 0.025))) +
      geom_line(aes(x = time/24, y = median, color=pk)) +
      geom_ribbon(aes(x=time/24, ymax = upper, ymin = lower, fill=pk), alpha = 0.3) +
      labs(x = "Time (days)", y = "Concentration (mg/L)", title = "Simulation Results") +
      xlim(c(0,100)) + theme_bw() + theme(legend.position = "top", legend.title = element_blank())
  })

  Cycle1 <- reactive({sim_data() %>%
    filter(time <= 24*input$ii_choice) %>%
    group_by(ID) %>%
    summarise(CMAXC1 = max(CMAX),
              CMINC1 = IPRED[n()],
              CAVGC1 = max(AUC)/24*input$ii_choice,
              CL = median(CL),
              Q = median(Q),
              V1 = median(V1),
              V2 = median(V2)) %>%
    ungroup() %>%
    mutate(QCMAXC1 = ntile(CMAXC1, 4),
           QCMINC1 = ntile(CMINC1, 4),
           QCAVGC1 = ntile(CAVGC1, 4))})

  TTEsim <- reactive({
    TTEmod <- mread("TTEmodel.cpp")
    TTEsim <- TTEmod %>%
      ev(amt = 100, ii = 24*input$ii_choice, addl = 50) %>%
      mrgsim_df(idata = Cycle1(), end = 24*input$sim_end, delta = 1) %>%
    group_by(ID) %>%
    mutate(RAND = runif(1)) %>%
    mutate(DV = ifelse(RAND >= ST, 1, 0),
           DV_PFS = ifelse(RAND >= ST_PFS, 1, 0))})

  TTEsim_null <- reactive({TTEsim() %>%
    group_by(ID) %>%
    filter(time == ifelse(sum(DV) >= 1, min(time[DV == 1]), max(time))) %>%
    left_join(Cycle1()) %>%
    ungroup() %>%
    mutate(QCMAX = ntile(CMAX, 4),
           QCMIN = ntile(CMIN, 4),
           QCAVG = ntile(CAVG, 4))})

  TTEsim_corr <- reactive({TTEsim() %>%
    group_by(ID) %>%
    select(-DV) %>%
    rename(DV = DV_PFS) %>%
    filter(time == ifelse(sum(DV) >= 1, min(time[DV == 1]), max(time))) %>%
    left_join(Cycle1()) %>%
    ungroup() %>%
    mutate(QCMAX = ntile(CMAX, 4),
           QCMIN = ntile(CMIN, 4),
           QCAVG = ntile(CAVG, 4))})

  pvalueqt <- reactive({
    fit.base0 <- coxph(Surv(time/24/30.475,DV)~1 , data=TTEsim_null())
    form <- as.formula(paste("Surv(time/24/30.475,DV)~",input$pkmetric,sep=""))
    sigmodel <- coxph(form , data=TTEsim_null())
    pvalueqt <-  signif(anova(sigmodel, fit.base0, test = "LRT")[2,4],3)
  })

  pvalueq <- reactive({
    fit.base0 <- coxph(Surv(time/24/30.475,DV)~1 , data=TTEsim_null())
    form <- as.formula(paste("Surv(time/24/30.475,DV)~",input$pkmetric,"C1",sep=""))
    sigmodel <- coxph(form , data=TTEsim_null())
    pvalueq <-  signif(anova(sigmodel, fit.base0, test = "LRT")[2,4],3)
  })

  pvalueqtcorr <- reactive({
    fit.base0 <- coxph(Surv(time/24/30.475,DV)~1 , data=TTEsim_corr())
    form <- as.formula(paste("Surv(time/24/30.475,DV)~",input$pkmetric,sep=""))
    sigmodel <- coxph(form , data=TTEsim_corr())
    pvalueqtcorr <-  signif(anova(sigmodel, fit.base0, test = "LRT")[2,4],3)
  })

  pvalueqcorr <- reactive({
    fit.base0 <- coxph(Surv(time/24/30.475,DV)~1 , data=TTEsim_corr())
    form <- as.formula(paste("Surv(time/24/30.475,DV)~",input$pkmetric,"C1",sep=""))
    sigmodel <- coxph(form , data=TTEsim_corr())
    pvalueqcorr <-  signif(anova(sigmodel, fit.base0, test = "LRT")[2,4],3)
  })

  output$km_plot_null <- renderPlot({
    surv_fit(Surv(time/24/30.475,DV)~1, data=TTEsim_null()) %>%
      ggsurvplot(xlab = "Time (Months)",
                 ylab = "Proportion without Event",
                 legend = "none",
                 break.x.by = 3,
                 palette = c("#39998e","#da674a","#7d3c98","#ffaa67","#90B0D9","#E5C76B")
      )
  })
    output$km_plot_nullqt <- renderPlot({
      form <- as.formula(paste("Surv(time/24/30.475,DV)~Q",input$pkmetric,sep=""))
      surv_fit(form, data=TTEsim_null()) %>%
        ggsurvplot(xlab = "Time (Months)",
                   ylab = "Proportion without Event",
                   title = paste("KM Curve by Quartiles of ",input$pkmetric," up to event time (p=",pvalueqt(),")", sep=""),
                   legend = "top",
                   break.x.by = 3,
                   palette = c("#39998e","#da674a","#7d3c98","#ffaa67","#90B0D9","#E5C76B")
        )
  })

    output$km_plot_nullq <- renderPlot({
      form <- as.formula(paste("Surv(time/24/30.475,DV)~Q",input$pkmetric,"C1",sep=""))
      surv_fit(form, data=TTEsim_null()) %>%
        ggsurvplot(xlab = "Time (Months)",
                   ylab = "Proportion without Event",
                   title = paste("KM Curve by Quartiles of Cycle 1 ",input$pkmetric," (p=",pvalueq(),")",sep=""),
                   legend = "top",
                   break.x.by = 3,
                   palette = c("#39998e","#da674a","#7d3c98","#ffaa67","#90B0D9","#E5C76B")
        )
    })

    output$km_plot_corr <- renderPlot({
      surv_fit(Surv(time/24/30.475,DV)~1, data=TTEsim_corr()) %>%
        ggsurvplot(xlab = "Time (Months)",
                   ylab = "Proportion without Event",
                   legend = "none",
                   break.x.by = 3,
                   palette = c("#39998e","#da674a","#7d3c98","#ffaa67","#90B0D9","#E5C76B")
        )
    })
    output$km_plot_corrqt <- renderPlot({
      form <- as.formula(paste("Surv(time/24/30.475,DV)~Q",input$pkmetric,sep=""))
      surv_fit(form, data=TTEsim_corr()) %>%
        ggsurvplot(xlab = "Time (Months)",
                   ylab = "Proportion with Event",
                   title = paste("KM Curve by Quartiles of ",input$pkmetric," up to event time (p=",pvalueqtcorr(),")", sep=""),
                   legend = "top",
                   break.x.by = 3,
                   palette = c("#39998e","#da674a","#7d3c98","#ffaa67","#90B0D9","#E5C76B")
        )
    })

    output$km_plot_corrq <- renderPlot({
      form <- as.formula(paste("Surv(time/24/30.475,DV)~Q",input$pkmetric,"C1",sep=""))
      surv_fit(form, data=TTEsim_corr()) %>%
        ggsurvplot(xlab = "Time (Months)",
                   ylab = "Proportion with Event",
                   title = paste("KM Curve by Quartiles of Cycle 1 ",input$pkmetric," (p=",pvalueqcorr(),")",sep=""),
                   legend = "top",
                   break.x.by = 3,
                   palette = c("#39998e","#da674a","#7d3c98","#ffaa67","#90B0D9","#E5C76B")
        )
    })
    newdata1 <- reactive({
      nameholder <- paste(input$pkmetric,"C1",sep="")

    newdata1 <- data.frame(ROWC = seq(1,10000,1))
    newdata1 <- newdata1 %>%
      mutate({{nameholder}} := seq(min(TTEsim_null()[,paste(input$pkmetric,"C1",sep="")]),max(TTEsim_null()[,paste(input$pkmetric,"C1",sep="")]),length.out = 10000))
    finalModel <- glm(as.formula(paste("DV ~ ",input$pkmetric,"C1",sep="")), family="binomial", data=TTEsim_null())
    preds <- predict(finalModel, newdata = newdata1, type = "link", se.fit = T)

    critval <- 1.96 ## approx 95% CI
    upr <- preds$fit + (critval * preds$se.fit)
    lwr <- preds$fit - (critval * preds$se.fit)
    fit <- preds$fit

    fit2 <- finalModel$family$linkinv(fit)
    upr2 <- finalModel$family$linkinv(upr)
    lwr2 <- finalModel$family$linkinv(lwr)

    newdata1$lwr <- lwr2
    newdata1$upr <- upr2
    newdata1$fit <- fit2
    newdata1 <- newdata1})

    LRpvalueq <- reactive({
      fit.base0 <- glm(as.formula(paste("DV ~ 1",sep="")), family="binomial", data=TTEsim_null())
      sigmodel <- glm(as.formula(paste("DV ~ ",input$pkmetric,"C1",sep="")), family="binomial", data=TTEsim_null())
      LRpvalueq <-  signif(anova(sigmodel, fit.base0, test = "LRT")[2,5],3)
    })

    summary <- reactive({TTEsim_null() %>%
      group_by(!!sym(paste("Q",input$pkmetric,"C1",sep=""))) %>%
      summarise(median = median(!!sym(paste(input$pkmetric,"C1",sep="")), na.rm = T),
                mean = mean(as.numeric(DV), na.rm = T),
                lower = mean + -qnorm(1-0.05/2)*sqrt((1/n())*mean*(1-mean)),
                upper = mean + qnorm(1-0.05/2)*sqrt((1/n())*mean*(1-mean)),
                Nlabel = paste(sum(DV==1),"/",n(), sep="")) %>%
      mutate(upper = ifelse(upper >= 1, 1, upper),
             lower = ifelse(lower <= 0, 0, lower))
    })
    output$lr_plot_nullq <- renderPlot({
      ggplot() + geom_jitter(data=TTEsim_null() %>% mutate(DV = ifelse(DV == 1,1.08,-0.08)), mapping=aes(x=!!sym(paste(input$pkmetric,"C1",sep="")),y=as.numeric(DV)), width = 0, height = 0.04, size = 1, alpha = 0.3) +
        geom_point(data=summary(), aes(x=median, y=mean), shape = 2) +
        geom_text(data=summary(), aes(x=median, y=upper+0.03, label = Nlabel), size = 3) +
        geom_errorbar(data=summary(), aes(x=median, ymin=lower, ymax=upper)) +
        geom_line(data=newdata1(), mapping=aes(x=!!sym(paste(input$pkmetric,"C1",sep="")), y=fit)) +
        geom_ribbon(data=newdata1(), aes(x=!!sym(paste(input$pkmetric,"C1",sep="")), ymin = lwr, ymax = upr), alpha = 0.2) +
        geom_text(data=summary() %>% slice(1), aes(x=median,y=0.8,label = paste("p = ",LRpvalueq(),sep=""))) +
        theme_bw() + xlab(paste("Cycle 1 ",input$pkmetric)) + ylab("Probability of Event") +
        scale_y_continuous(breaks = seq(0,1,0.2))

    })

    newdata1t <- reactive({
      nameholder <- paste(input$pkmetric,sep="")

      newdata1 <- data.frame(ROWC = seq(1,10000,1))
      newdata1 <- newdata1 %>%
        mutate({{nameholder}} := seq(min(TTEsim_null()[,paste(input$pkmetric,sep="")]),max(TTEsim_null()[,paste(input$pkmetric,sep="")]),length.out = 10000))
      finalModel <- glm(as.formula(paste("DV ~ ",input$pkmetric,sep="")), family="binomial", data=TTEsim_null())
      preds <- predict(finalModel, newdata = newdata1, type = "link", se.fit = T)

      critval <- 1.96 ## approx 95% CI
      upr <- preds$fit + (critval * preds$se.fit)
      lwr <- preds$fit - (critval * preds$se.fit)
      fit <- preds$fit

      fit2 <- finalModel$family$linkinv(fit)
      upr2 <- finalModel$family$linkinv(upr)
      lwr2 <- finalModel$family$linkinv(lwr)

      newdata1$lwr <- lwr2
      newdata1$upr <- upr2
      newdata1$fit <- fit2
      newdata1t <- newdata1})

    LRpvalueqt <- reactive({
      fit.base0 <- glm(as.formula(paste("DV ~ 1",sep="")), family="binomial", data=TTEsim_null())
      sigmodel <- glm(as.formula(paste("DV ~ ",input$pkmetric,sep="")), family="binomial", data=TTEsim_null())
      LRpvalueqt <-  signif(anova(sigmodel, fit.base0, test = "LRT")[2,5],3)
    })

    summaryt <- reactive({TTEsim_null() %>%
        group_by(!!sym(paste("Q",input$pkmetric,sep=""))) %>%
        summarise(median = median(!!sym(paste(input$pkmetric,sep="")), na.rm = T),
                  mean = mean(as.numeric(DV), na.rm = T),
                  lower = mean + -qnorm(1-0.05/2)*sqrt((1/n())*mean*(1-mean)),
                  upper = mean + qnorm(1-0.05/2)*sqrt((1/n())*mean*(1-mean)),
                  Nlabel = paste(sum(DV==1),"/",n(), sep="")) %>%
        mutate(upper = ifelse(upper >= 1, 1, upper),
               lower = ifelse(lower <= 0, 0, lower))
    })
    output$lr_plot_nullqt <- renderPlot({
      ggplot() + geom_jitter(data=TTEsim_null() %>% mutate(DV = ifelse(DV == 1,1.08,-0.08)), mapping=aes(x=!!sym(paste(input$pkmetric,sep="")),y=as.numeric(DV)), width = 0, height = 0.04, size = 1, alpha = 0.3) +
        geom_point(data=summaryt(), aes(x=median, y=mean), shape = 2) +
        geom_text(data=summaryt(), aes(x=median, y=upper+0.03, label = Nlabel), size = 3) +
        geom_errorbar(data=summaryt(), aes(x=median, ymin=lower, ymax=upper)) +
        geom_line(data=newdata1t(), mapping=aes(x=!!sym(paste(input$pkmetric,sep="")), y=fit)) +
        geom_ribbon(data=newdata1t(), aes(x=!!sym(paste(input$pkmetric,sep="")), ymin = lwr, ymax = upr), alpha = 0.2) +
        geom_text(data=summaryt() %>% slice(1), aes(x=median,y=0.8,label = paste("p = ",LRpvalueqt(),sep=""))) +
        theme_bw() + xlab(paste(input$pkmetric," up to event time",sep="")) + ylab("Probability of Event") +
        scale_y_continuous(breaks = seq(0,1,0.2))

    })

    newdata1corr <- reactive({
      nameholder <- paste(input$pkmetric,"C1",sep="")

      newdata1 <- data.frame(ROWC = seq(1,10000,1))
      newdata1 <- newdata1 %>%
        mutate({{nameholder}} := seq(min(TTEsim_corr()[,paste(input$pkmetric,"C1",sep="")]),max(TTEsim_corr()[,paste(input$pkmetric,"C1",sep="")]),length.out = 10000))
      finalModel <- glm(as.formula(paste("DV ~ ",input$pkmetric,"C1",sep="")), family="binomial", data=TTEsim_corr())
      preds <- predict(finalModel, newdata = newdata1, type = "link", se.fit = T)

      critval <- 1.96 ## approx 95% CI
      upr <- preds$fit + (critval * preds$se.fit)
      lwr <- preds$fit - (critval * preds$se.fit)
      fit <- preds$fit

      fit2 <- finalModel$family$linkinv(fit)
      upr2 <- finalModel$family$linkinv(upr)
      lwr2 <- finalModel$family$linkinv(lwr)

      newdata1$lwr <- lwr2
      newdata1$upr <- upr2
      newdata1$fit <- fit2
      newdata1corr <- newdata1})

    LRpvalueqcorr <- reactive({
      fit.base0 <- glm(as.formula(paste("DV ~ 1",sep="")), family="binomial", data=TTEsim_corr())
      sigmodel <- glm(as.formula(paste("DV ~ ",input$pkmetric,"C1",sep="")), family="binomial", data=TTEsim_corr())
      LRpvalueqcorr <-  signif(anova(sigmodel, fit.base0, test = "LRT")[2,5],3)
    })

    summarycorr <- reactive({TTEsim_corr() %>%
        group_by(!!sym(paste("Q",input$pkmetric,"C1",sep=""))) %>%
        summarise(median = median(!!sym(paste(input$pkmetric,"C1",sep="")), na.rm = T),
                  mean = mean(as.numeric(DV), na.rm = T),
                  lower = mean + -qnorm(1-0.05/2)*sqrt((1/n())*mean*(1-mean)),
                  upper = mean + qnorm(1-0.05/2)*sqrt((1/n())*mean*(1-mean)),
                  Nlabel = paste(sum(DV==1),"/",n(), sep="")) %>%
        mutate(upper = ifelse(upper >= 1, 1, upper),
               lower = ifelse(lower <= 0, 0, lower))
    })
    output$lr_plot_corrq <- renderPlot({
      ggplot() + geom_jitter(data=TTEsim_corr() %>% mutate(DV = ifelse(DV == 1,1.08,-0.08)), mapping=aes(x=!!sym(paste(input$pkmetric,"C1",sep="")),y=as.numeric(DV)), width = 0, height = 0.04, size = 1, alpha = 0.3) +
        geom_point(data=summarycorr(), aes(x=median, y=mean), shape = 2) +
        geom_text(data=summarycorr(), aes(x=median, y=upper+0.03, label = Nlabel), size = 3) +
        geom_errorbar(data=summarycorr(), aes(x=median, ymin=lower, ymax=upper)) +
        geom_line(data=newdata1corr(), mapping=aes(x=!!sym(paste(input$pkmetric,"C1",sep="")), y=fit)) +
        geom_ribbon(data=newdata1corr(), aes(x=!!sym(paste(input$pkmetric,"C1",sep="")), ymin = lwr, ymax = upr), alpha = 0.2) +
        geom_text(data=summarycorr() %>% slice(1), aes(x=median,y=0.8,label = paste("p = ",LRpvalueqcorr(),sep=""))) +
        theme_bw() + xlab(paste("Cycle 1 ",input$pkmetric)) + ylab("Probability of Event") +
        scale_y_continuous(breaks = seq(0,1,0.2))

    })

    newdata1tcorr <- reactive({
      nameholder <- paste(input$pkmetric,sep="")

      newdata1 <- data.frame(ROWC = seq(1,10000,1))
      newdata1 <- newdata1 %>%
        mutate({{nameholder}} := seq(min(TTEsim_corr()[,paste(input$pkmetric,sep="")]),max(TTEsim_corr()[,paste(input$pkmetric,sep="")]),length.out = 10000))
      finalModel <- glm(as.formula(paste("DV ~ ",input$pkmetric,sep="")), family="binomial", data=TTEsim_corr())
      preds <- predict(finalModel, newdata = newdata1, type = "link", se.fit = T)

      critval <- 1.96 ## approx 95% CI
      upr <- preds$fit + (critval * preds$se.fit)
      lwr <- preds$fit - (critval * preds$se.fit)
      fit <- preds$fit

      fit2 <- finalModel$family$linkinv(fit)
      upr2 <- finalModel$family$linkinv(upr)
      lwr2 <- finalModel$family$linkinv(lwr)

      newdata1$lwr <- lwr2
      newdata1$upr <- upr2
      newdata1$fit <- fit2
      newdata1tcorr <- newdata1})

    LRpvalueqtcorr <- reactive({
      fit.base0 <- glm(as.formula(paste("DV ~ 1",sep="")), family="binomial", data=TTEsim_corr())
      sigmodel <- glm(as.formula(paste("DV ~ ",input$pkmetric,sep="")), family="binomial", data=TTEsim_corr())
      LRpvalueqtcorr <-  signif(anova(sigmodel, fit.base0, test = "LRT")[2,5],3)
    })

    summarytcorr <- reactive({TTEsim_corr() %>%
        group_by(!!sym(paste("Q",input$pkmetric,sep=""))) %>%
        summarise(median = median(!!sym(paste(input$pkmetric,sep="")), na.rm = T),
                  mean = mean(as.numeric(DV), na.rm = T),
                  lower = mean + -qnorm(1-0.05/2)*sqrt((1/n())*mean*(1-mean)),
                  upper = mean + qnorm(1-0.05/2)*sqrt((1/n())*mean*(1-mean)),
                  Nlabel = paste(sum(DV==1),"/",n(), sep="")) %>%
        mutate(upper = ifelse(upper >= 1, 1, upper),
               lower = ifelse(lower <= 0, 0, lower))
    })
    output$lr_plot_corrqt <- renderPlot({
      ggplot() + geom_jitter(data=TTEsim_corr() %>% mutate(DV = ifelse(DV == 1,1.08,-0.08)), mapping=aes(x=!!sym(paste(input$pkmetric,sep="")),y=as.numeric(DV)), width = 0, height = 0.04, size = 1, alpha = 0.3) +
        geom_point(data=summarytcorr(), aes(x=median, y=mean), shape = 2) +
        geom_text(data=summarytcorr(), aes(x=median, y=upper+0.03, label = Nlabel), size = 3) +
        geom_errorbar(data=summarytcorr(), aes(x=median, ymin=lower, ymax=upper)) +
        geom_line(data=newdata1tcorr(), mapping=aes(x=!!sym(paste(input$pkmetric,sep="")), y=fit)) +
        geom_ribbon(data=newdata1tcorr(), aes(x=!!sym(paste(input$pkmetric,sep="")), ymin = lwr, ymax = upr), alpha = 0.2) +
        geom_text(data=summarytcorr() %>% slice(1), aes(x=median,y=0.8,label = paste("p = ",LRpvalueqtcorr(),sep=""))) +
        theme_bw() + xlab(paste(input$pkmetric," up to event time",sep="")) + ylab("Probability of Event") +
        scale_y_continuous(breaks = seq(0,1,0.2))

    })


    # Add the sim_data() reactive expression
    output$sim_table <- renderTable({
      TTEsim_null()
      # Add the sim_data() reactive expression
    })
}

shinyApp(ui = ui, server = server)

