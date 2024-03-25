#Run simulation for AERA 2024
#comparing algorithms to select number of classes

#install.packages(c("GA", "VarSelLCM", "LCAvarsel", "poLCA" ))

#load packages
library(GA) #for genetic algorithm
library(VarSelLCM)
library(LCAvarsel)
library(poLCA)


start.time = proc.time()

#set parameters of informative items
#these parameters were chosen through trial and error so that
#with a population of 100000, entropy was 0.8 (calculated with Mplus)
all.tau= list(tau2class = rbind(class1= c(-1,-1.5,-1,-0.5,0,-0.8,-1.2,-1.5,-0.7,-0.4), #high ability
                                class2 = c(0.9, 1.2, 2,1,1.5,1,0.4,1.8,0.7,1.3)), #low ability
                                
                
              tau3class = rbind(class1= c(-2,-1.7,-1.5,-1.3,-1, -0.7, -0.3, -1, -0.3, -0.5), #high ability
                            class2 = c(1, -2, 1, -1, 2, -1, 1, -1.3, 1.5, -2.5), #mixed ability
                            class3 = c(2,1.3, 1.5, 2,1.5,2,1.8, 1.7, 2, 1.5))) #low ability


#conditions
#entropy = 0.8
#number of informative items = 5
#number of redundant items: 5,10
#number of non-informative items: 20,100
#sample size = 20000, 40000
all.results = data.frame() #storage
iters = 2 #total number of iterations

for (classes in 2:3) {

  sample = 2000


 if (classes==2) {tau= all.tau[[1]]} else {tau=all.tau[[2]]} 

#20 redundant items are two copies of tau
redundant.par = cbind(tau,tau)

#20 non.informative (irrelevant) items have randomly generated parameters that are equal across classes
non.inf.par = matrix(rep(runif(20, min=min(tau), 
                               max=max(tau)),each=nrow(tau)),nrow=nrow(tau),
                     ncol=20)
#final parameters
par = cbind(tau,redundant.par,non.inf.par)


#assign class sizes
if (classes == 2 ) {sizes =  c(0.4*sample, 0.6*sample)}
if (classes == 3 ) {sizes = c(0.3*sample, 0.3*sample, 0.4*sample)}
#===================================================
#SIMULATE DATA

#convert logit to probability
prob <- 1/(1+exp(par))

#convert probabilities to a list of matrix required by poLCA.simdata
all.probs = list()
for (c in 1:ncol(prob)) {
  prob.matrix = cbind(prob[,c],1-(prob[,c]))
  all.probs[[c]] = prob.matrix
}


#loop through iterations
for (i in 1:iters) {

  #simulate data
  all.data = poLCA.simdata(N=sum(sizes), probs=all.probs, P=(sizes/sum(sizes)))
  all.data = all.data$dat
  

#fit LCA 
#require(poLCA)
#LCA.model1 = poLCA(cbind(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10)~1,Data,nclass=2,calc.se = F)

#loop through redundant and irrelevant
for (redundant in 1:3) {
  for (irrelevant in 1:3) {

n.redundant = list(NULL,11:20,11:30)[[redundant]]
n.irrelevant = list(NULL,31:40,31:50)[[irrelevant]]

    
#select subsets of the dataset with appropriate redundant and irrelevant
  data = all.data[,c(1:10,n.redundant, n.irrelevant)]
  for (var in 1:ncol(data)) {data[,var]=as.factor(data[,var])}

#results for VarSelLCM (Marback and Sedki, 2017)
results1 = VarSelCluster(data,gvals=1:5,vbleSelec = T,crit.varsel = "BIC")

#results for LCAvarsel (Fop, Smart, and Murphy) with genetic algorithm
results2 = LCAvarsel(data, G=1:5, search="ga", verbose=F)

#results for LCAvarsel (Fop, Smart, and Murphy) with backward stepwise search algorithm
results3 = LCAvarsel(data, G=1:5, search="backward",swap=T, verbose=F)

#results for LCAvarsel (Fop, Smart, and Murphy) with forward stepwise search algorithm
results4 = LCAvarsel(data, G=1:5, search="forward", swap=T, verbose=F)
#---------------------------------------

#count the number of classes extracted
classes1 = results1@model@g
classes2 = results2$model$G
classes3 = results3$model$G
classes4 = results4$model$G

#--------------------------------------
# #count the number of relevant items selected (which includes redundant)
relevant1 = length(intersect(paste("Y",1:30,sep=""),
                      results1@model@names.relevant))
relevant2 = length(intersect(paste("Y",1:30,sep=""),
                             results2$variables))
relevant3 = length(intersect(paste("Y",1:30,sep=""),
                             results3$variables))
relevant4 = length(intersect(paste("Y",1:30,sep=""),
                             results4$variables))
#-----------------------------------------------------------
#count the number of irrelevant items selected
irrelevant1 = length(intersect(paste("Y",31:50,sep=""),
                             results1@model@names.relevant))
irrelevant2 = length(intersect(paste("Y",31:50,sep=""),
                             results2$variables))
irrelevant3 = length(intersect(paste("Y",31:50,sep=""),
                             results3$variables))
irrelevant4 = length(intersect(paste("Y",31:50,sep=""),
                             results4$variables))

#------------------------------------------
#count redundant items selected
redundant.table = matrix(paste("Y",1:30,sep=""),10,3)
redundant1 =c()
redundant2 = c()
redundant3 = c()
redundant4 = c()

for (r in 1:nrow(redundant.table)) {
  #count the number in each row minus 1 (any more than 1 are redundant)
 redundant1 = c(redundant1,length(intersect(redundant.table[r,],
                                               results1@model@names.relevant))) 
 redundant2 = c(redundant2,length(intersect(redundant.table[r,],
                                                         results2$variables))) 
 redundant3 = c(redundant3,length(intersect(redundant.table[r,],
                                            results3$variables)))
 redundant4 = c(redundant4,length(intersect(redundant.table[r,],
                                            results4$variables)))
 }
redundant1 =sum(redundant1>1)
redundant2 =sum(redundant2>1)
redundant3 =sum(redundant3>1)
redundant4 =sum(redundant4>1)

#------------------------------ 
# #put together results
condition.results = data.frame(i,classes,redundant, irrelevant,
                          classes1,classes2,
                         relevant1,relevant2, 
                         irrelevant1,irrelevant2,
                         redundant1, redundant2)

all.results = rbind(all.results, condition.results)
# 
# 
# #save current set of results
# write.table(selected.results, file="simulation_results.csv", row.names = F, 
#             col.names=F, append=T, sep=",")



#close loops
} } 



# #label results
# selected.results = data.frame(current.condition, selected.results)
# all.results = rbind(all.results, selected.results)
                                       

} #close loop through iterations




} #close loop through classes

end.time=proc.time()
end.time - start.time

#write all results
write.table(all.results, "all_results.csv", row.names = F, quote=F, sep=",")

