animals = c("706", "716", "721", "724R", "727R")
treatments = c("CN", "BCG", "MAP", "BOVIS")
timepoints = c("2", "6", "24")

values = c(1:60)
f_timepoints = as.factor(c(rep("2",20), rep("6",20), rep("24",20)))
f_treatments = as.factor(rep(c(rep("CN",5),rep("BCG",5),rep("MAP",5),rep("BOVIS",5)),3))
f_animals = as.factor(rep(rep(c("706", "716", "721", "724R", "727R"),4),3))

glm.add = glm(values ~ f_animals + f_timepoints + f_treatments, family=gaussian())
anova(glm.add)
summary(glm.add)

glm.multi = glm(values~f_animals * f_timepoints * f_treatments, family=gaussian())anova(glm.fake)
anova(glm.multi) # not enough degrees of freedom (only one replicate for interaction effects)
summary(glm.multi) # NAs for Std. Error because only one value per category

# Missing f_animals706 from the Coefficients, I think I forgot a constant in the model
glm.add.constant = glm(values ~ 0 + f_animals + f_timepoints + f_treatments, family=gaussian())
anova(glm.add.constant)
summary(glm.add.constant) # Missing 2HR time and intercept estimate

# Not sure whether the constant should be 1 or 0
glm.add.constant.one = glm(values ~ 1 + f_animals + f_timepoints + f_treatments, family=gaussian())
anova(glm.add.constant.one)
summary(glm.add.constant.one) # Now there is an intercept estimate again