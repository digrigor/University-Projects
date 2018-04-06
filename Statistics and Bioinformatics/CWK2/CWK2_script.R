#Name: Dionysios Grigoriadis
#ID: 170484367
#BIO782P Statistics and Bioinformatics Assignment 2017
#CWK 2
#Part 1

stimulation_potential = seq(-100, 40, by=10)
response_voltage_control = c(-70.00000,-70.00000,-70.00000,-70.00000,-70.00000,35.79871,38.75082,39.30423,39.65390,35.55835,35.66762,36.55754,39.51027,35.65274,39.17013)
response_voltage_concN_1 = c(-70.00000,-70.00000,-70.00000,-70.00000,37.88144,35.64345,36.54319,36.18773,36.67542,38.56634,38.53102,37.73421,38.69826,38.67744,37.21896)
response_voltage_concN_2 = c(-70.00000,-70.00000,-70.00000,37.16110,39.28471,36.57506,38.03037,39.37853,39.92131,39.01116,39.46393,36.16073,37.70418,35.88598,38.53137)
response_voltage_concN_3 = c(-70.00000,-70.00000,36.25729,35.70035,35.61478,36.81245,37.01299,39.34139,36.59520,35.00749,37.94440,36.10611,37.15617,35.29182,35.32412)

plot(stimulation_potential,stimulation_potential,type='n',xlab='Stimulation, mV',ylab='Excitation, mV')
lines(stimulation_potential,response_voltage_control)
lines(stimulation_potential,response_voltage_concN_1,col="green")
lines(stimulation_potential,response_voltage_concN_2,col="orange")
lines(stimulation_potential,response_voltage_concN_3,col="red")

legend(0,20,"control",fill="black")
legend(0,00,"concentration 1",fill="green")
legend(0,-20,"concentration 2",fill="orange")
legend(0,-40,"concentration 3",fill="red")


predict = function(activation_threshold,stimulation_voltage){
  activation = rep(-70,length(stimulation_voltage))
  activation[stimulation_voltage>activation_threshold] = 40
  return(activation)
}
calculate_errors = function(predicted, observed){
  total_errors = sum((observed - predicted)^2)
  return(total_errors)
}
fit_threshold = function(input_values_stimulation,input_values_response,threshold){
  predicted_values = predict(threshold,input_values_stimulation)
  fit_errors = calculate_errors(predicted_values,input_values_response)
  return(fit_errors)
}

#Question 6:
#The existing optimisation method has some restrictions:

fitted_threshold = runif(1,-100,40)
error = fit_threshold(stimulation_potential,response_voltage_control,fitted_threshold)

for(i in 1:20){ #Low numbers of times that the optimisation loop runs, leading to
  #lower possibility of finding the best possible prediction.
  new_threshold = runif(1,-100,40) #The proposal mechanisms keeps proposing 
  #a random threshold each time, even if the previous guess was close to the true value.
  new_error     = fit_threshold(stimulation_potential,response_voltage_control,new_threshold)
  if(new_error < error){
    fitted_threshold = new_threshold
    error = new_error
  }
}
fitted_threshold
error

#Proposing a better optimisation method:

new_threshold = runif(1,-100,40)
for(i in 1:1000){ #The number of times that the optimisation loop runs is increased 
  new_error = fit_threshold(stimulation_potential,response_voltage_control,new_threshold)
  if(new_error < error){
    fitted_threshold = new_threshold
    error = new_error
    new_threshold = runif(1, fitted_threshold-40, fitted_threshold+40) #If the optimiser finds a better 
    #threshold, the next threshold will be picked within a range around the previous guess.
  }
  else if(new_error >= error){
    new_threshold = runif(1,-100,40) #The new threshold produces a worse prediction
    #compared to the previous one: A new random threshold is picked for the next 
    #optimisation round.
  }
  }
fitted_threshold
error

write.csv(data.frame('fitted value'=fitted_threshold,'errors SS'=error),file = 'excitation.csv',row.names = F)

#Question 7:
#The code was modified in order to fit the data given under the control experiment
#and the three experimental nicotine treatments. 

response_list=list(response_voltage_control, response_voltage_concN_1, response_voltage_concN_2, response_voltage_concN_3)
errors=c() #List of all errors produced by the optimisations
fitted_thresholds=c() #List of all thresholds produced by the optimisations

for(count in 1:10){#The optimization loop will run 10 times for each experiment (control,
  #concentration 1, concnetration 2, concentration 3)
  
  for(z in response_list){#For an experiment each time.
    
    fitted_threshold = runif(1,-100,40)
    error = fit_threshold(stimulation_potential,z,fitted_threshold) #A model is fitted
    #under a z experimental treatment (control,concentration 1, concnetration 2, 
    #concentration 3). The goodness of fit is calculated. The optimisation method
    #described previously is running.
    new_threshold = runif(1,-100,40)
    for(i in 1:100){
      new_error = fit_threshold(stimulation_potential,z,new_threshold) 
      
      if(new_error < error){
        fitted_threshold = new_threshold
        error = new_error
        new_threshold = rnorm(1, fitted_threshold, 30)
      }
      else if(new_error >= error){
        new_threshold = runif(1,-100,40)
      }
    }
    errors[(length(errors) + 1)] <- error #List contatining all the error values
    #for each optimisation that has ran for all the 4 experiments.
    fitted_thresholds[(length(fitted_thresholds) + 1)] <- fitted_threshold
    #List contatining all the fitted thresholds for each optimisation that has 
    #ran for all the 4 experiments. 
    #In these lists, the values follow the repetitive pattern (control, concN 1, concN 2,
    #concN 3 for 10 times).
  }
}
#Comparing the goodness of fit between the models: 10 error values are produced by the 
#optimisation method for the control experiment and each experimental Nicotine treatment.
#The mean error value for each experiment is calculated and stored into the final
#compare_errors matrix.
errors_table=data.frame(matrix(errors,nrow=10 ,ncol=4, byrow = TRUE))
colnames(errors_table)=c("Control Err.", "Concentration1 Err.", "Concentration2 Err.", "Concentration3 Err.")
#Building the final comperative table
compare_errors=matrix(colMeans(errors_table))
rownames(compare_errors)=c("Control", "Concentration 1", "Concentration 2", "Concentration 3")
print(compare_errors)

#Comparing the fitted activation threshold between the control and the experimental 
#treatments: 10 fitted activation threshold values are produced by the optimisation 
#method for the control experiment and each experimental Nicotine treatment.
fit_thresholds_table=data.frame(matrix(fitted_thresholds,nrow=10 ,ncol=4, byrow = TRUE))
colnames(fit_thresholds_table)=c("Control", "Concentration 1", "Concentration 2", "Concentration 3")
attach(fit_thresholds_table)
#The threshold values of the fitted models of the control values are
#compared with the respective values of each experimental treatment with nicotine.
test1=t.test(`Control`,`Concentration 1`)
test2=t.test(`Control`, `Concentration 2`)
test3=t.test(`Control`, `Concentration 3`)

par(mar=c(5, 6, 4, 2) + 0.1) #Extends the plot margins
boxplot(fit_thresholds_table, ylab="Predicted Activation Threshold\n(mV)")

#Building the final comperative table of the fitted activation thresholds
#between experimental treatments and control experiment
compare_thresholds=matrix(data = c(test1$p.value,test2$p.value,test3$p.value,
                                   mean(`Concentration 1`)/mean(Control),
                                   mean(`Concentration 2`)/mean(Control),
                                   mean(`Concentration 3`)/mean(Control)), nrow = 3, ncol = 2)
row.names(compare_thresholds)=c("Control/concN1","Control/concN2","Control/concN3")
colnames(compare_thresholds)=c("P-value","Fold Change")
print(compare_thresholds)

#The comparative table of the fitted activation thresholds between experimental 
#treatments and control experiment suggests that each nicotine treatment lowers 
#the activation threshold (see fold change values and pvalues in the table). The 
#null hypothesis in this case is that the activation threshold caused by each nicotine
#treatment is similar to the activation threshold observed on the control axons, thus
#the nicotine concentration does not affect the activation threshold of the axons. 
