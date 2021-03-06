# Import data
gene_data <- read.delim("test-data-skoog.txt")
# ?read.delim()
head(gene_data)
# Data are representative of taxa as rows and columns showing genes (values equal copies of those genes for each taxa). These 4 genes are required for the full pathway. We will show some examples of how to subset, filter, and perform basic calculations.

# Colnames and rownames, "index" what appears in the columns and rows. It is like the bar A-Z at the top and the numbers on the left side bar of an Excel table.
colnames(gene_data) # reports all of my column names! 

# Let's change some of these because there is a typo
colnames(gene_data)[6] <- "GENE_D" #Changes the 6th column name to "GENE_D"
# Change columns headers for columns 4-6
# colnames(gene_data)[4:6] <- c("GENE_C", "GENE_B", "GENE_D")

colnames(gene_data) # Check output

# Example of how to change the order (2 ways) - Which way would be more reproducible?
colnames(gene_data)
# Option 1
gene_data_reorder <- gene_data[, c("Taxon", "MAG","GENE_A", "GENE_B", "GENE_C", "GENE_D")]

# Option 2
tmp <- gene_data[c(1:3,5,4,6)] #Combine specific columns "c()"



# How about rows - these are just random numbers! R assigns increasing numbers automatically as row names, unless you specify. 
row.names(gene_data) # So when would you need to specify row.names???
# Recall the data I'm looking at are integers and factors. One row states the full taxon name and another row tells me the MAG ID number.
class(gene_data)
str(gene_data)

# In what scenario is it good to add in row.names? It is the same as adding column headers, but this is with respect to rows. Keeping with that logic, row names need to be unique.
tmp <- gene_data
row.names(tmp) <- tmp$Taxon
row.names(tmp) <- tmp$MAG
# head(tmp)
# View(tmp)
# row.names(tmp) <- paste(tmp$Taxon, tmp$MAG, sep=" ")


# Let's work on modifying the data frame and do some math with dplyr.
# library(reshape2)
library(tidyverse) # plyr # dplyr
# head(gene_data[1:2,])

# Common base R commands to check out help menus for:
# ?subset()
# ?apply() # Also tapply(), sapply(), and lapply()
# ?aggregate()
# ?split()


# Let's subset our data with subset()
# I only want MAGs (rows or entries) where Gene A is great that zero
head(gene_data[1:3,])
tmp1 <- subset(gene_data, !(GENE_A == 0 ))
tmp1 <- subset(gene_data, GENE_A > 0)

## Other approaches
# tmp1 <- subset(gene_data, (GENE_C > 0 & GENE_D > 0))
# tmp1 <- subset(gene_data, MAG %in% "MAG_2")

# Let's repeat this with tidyverse
# ?filter()
tmp2 <- filter(gene_data, GENE_A > 0)

# What if I want to filter a set of MAGs of interest?
mags_i_want <- c("MAG_4", "MAG_12", "MAG_24", "MAG_18")

tmp3 <- filter(gene_data, MAG %in% mags_i_want)

# How would we combine the statements above?


# Pipe operator - the key to tidy syntax!
gene_data %>% filter(MAG %in% mags_i_want)
head(gene_data)
gene_data %>% head()


# Filter to the mags I want, but only show the taxonomic identities. Also set it equal to a new data frame.
gene_data %>% 
  filter(MAG %in% mags_i_want) %>%
  select(Taxon, MAG)



# Repeat in base R
subset(gene_data, MAG %in% mags_i_want)[,c("Taxon", "MAG")]


# This data is in wide format. Convert to long format - this is more versatile in R.
## Compare long vs. wide format
# ?pivot_longer()
# ?pivot_wider()
colnames(gene_data)
gene_data_long <- gene_data %>%
  pivot_longer(GENE_A:GENE_D, names_to = "GENE", values_to = "COUNT")


# View(gene_data)
# View(gene_data_long)


# Split Taxon column
gene_data_long_cols <- gene_data_long %>%
  separate(Taxon, c("phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = FALSE)
head(gene_data_long_cols)


#Lets combine the two operations from above.
gene_data_long <- gene_data %>%
  pivot_longer(GENE_A:GENE_D, names_to = "GENE", values_to = "COUNT") %>%
  separate(Taxon, c("phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = FALSE) %>%
  data.frame
# head(gene_data_long)

# Mean at Class level - across whole dataset (add Min and Max value to demonstrate syntax)
# ?summarise()
mean_class <- gene_data_long %>%
  group_by(class, GENE) %>%
  summarise(MEAN = mean(COUNT)) %>%
  data.frame
# head(mean_class)


# Relative abundance of genes by Taxon
# ?mutate()
relabun_all <- gene_data_long %>%
  group_by(class, GENE) %>%
  mutate(RELABUN = (COUNT/sum(COUNT))) %>%
  data.frame
head(relabun_all)

## Check relative abundance calculation
# alpha_test <- filter(relabun_all, class %in% "c__Alphaproteobacteria")
# View(alpha_test)
# sum(alpha_test$RELABUN)


# Add additional calculations to the above summarise and mutate functions.
# ?n()
# ?n_distinct()
# ?median()
# ?max()
# ?min()
# ?var()
# ?sd()


# IF we have time to cover:

# Average across replicates
## Let's simulate replicates that we need to average across
head(gene_data_long)
# Add in a simulated column that represents replicates
gene_data_long$REPLICATE <- sample(c("rep-1", "rep-2"), replace = TRUE)
head(gene_data_long)

gene_avg <- gene_data_long %>%
  group_by(Taxon, MAG, variable) %>%
  summarise(AVG_REPS = mean(value)) %>%
  data.frame
head(gene_avg)  
# dim(gene_avg); dim(gene_data_long)


# Export Table for use
head(mean_class)
?write.table()



# Addressing questions from R slack Channel

# Change to binary
# gene_data_long$PresenceAbsence<-ifelse(gene_data_long$value > 0, 1, 0)
#
# Review from last week - explaination:
# iris[order(iris$Species, -iris$Petal.Length),]
#
# With "[ ]" we use a comma to designate if rows or columns are being considered within the brackets. Without a comma, it defaults to considering columns.
# By default it is: [rows, columns]
#
# Let's show an example using "head()"
head(iris) #Prints the top section of your data frame (6 lines total)
# Let's play around with adding brackets
head(iris[1,]) # prints the 1st row of the data frame
head(iris[1:3,]) #prints the 1st 3 rows of the data frame
## ^That's what the colon is for! it selects colums 1 through 3 (inclusive)
## If we move the comma
head(iris[,1:3]) #Now we print the 1st 3 columns only
head(iris[1:3]) # This also prints the 1st 3 columns - remember that R defaults to columns when no comma is present.
#
# Print the 1st 3 rows of the first 2 columns
head(iris[1:3,1:2])
#
head(iris)
# Let's break this function down:
# "order(iris$Species, -iris$Petal.Length)" will order the iris data frame by Species first (ascending) and THEN by Petal.length (descending because of the - sign). This is followed by a comma because we are ordering the rows.
iris[order(iris$Species, -iris$Petal.Length), ]
#
# Let's try adding something! 
iris[order(iris$Species, -iris$Petal.Length), 1:3]