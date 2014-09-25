################################################################################################################
# Program creates a list of parameter that are significant for producing a particular outcome                  #
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Braun lab group, Biology Department, University of Florida}      #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
###############################################################################################################



data_matrix <- read.csv('file.csv')
data_matrix <- data_matrix[,colSums(is.na(data_matrix)) != nrow(data_matrix)]

PcolNum <- ncol(data_matrix) - 1

initColNames <- colnames(data_matrix)[-1]

colnames(data_matrix) <- c("Outcome", paste("V", 1:PcolNum, sep=""))

total_parameter <- (ncol(data_matrix)-1)

outcome <- data_matrix$Outcome

D <- numeric(length(data_matrix)-1)

for (j in 1:length(data_matrix)-1) {
    D[j] <- paste("V", j, sep="")
}

total_data = ''
reduced_data=''

l=1


while(l <= total_parameter){
	total_data <- union(total_data, list(eval(parse(text=paste("data_matrix$", D[l], sep="")))))
    l=l+1
    }

total_data <- total_data[-1]

valAdd <- paste(total_data, collapse="+")
formula <- paste("outcome ~ ", valAdd, sep = "")
full <- lm(formula)

omit_parameter = 1

pvalList <- new.env()

while (omit_parameter <= total_parameter) {
	reduced_data <- total_data[-omit_parameter]
	reduced_valAdd <- paste(reduced_data, collapse="+")
	reduced_formula <- paste("outcome ~ ", reduced_valAdd, sep = "")
	reduced <- lm(reduced_formula)
	output <- anova(reduced, full)
        pvalList[[paste("V", omit_parameter, sep="")]] = output[[6]][2]
        omit_parameter <- omit_parameter + 1
}


counter = 1

while (counter <= length(pvalList)) {
    if (pvalList[[paste("V", counter, sep="")]] <= 0.05) {
        capture.output(print(paste0(initColNames[counter], 'is important \n')), file = 'out.txt', append = TRUE)
        
    }
    counter = counter + 1
}
