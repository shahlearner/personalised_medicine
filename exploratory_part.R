rm(list=ls())
setwd("F:/shah/personalized medicines/dirctry")

#Load libraries  

library('readr') # data input
library('tibble') # data wrangling
library('dplyr') # data manipulation
library('tidyr') # data wrangling
library('corrplot') # visualisation
library('ggfortify') # visualisation
library('grid') # visualisation
library('gridExtra') # visualisation
library('ggraph') # visualisation
library('ggplot2') # visualization
library('stringr') # string manipulation
library('forcats') # factor manipulation
library('tidytext') # text mining
library('SnowballC') # text analysis
library('wordcloud') # test visualisation

#load data files
#Reading in the variants data tables
#train variants
train_data = read.csv("train_variants_csv.csv",header = T)
test_data = read.csv("test_variants_csv.csv",header = T)


#loading text files
#training_text
train_text_org = tibble(text = read_lines('training_text', skip = 1))#reading text file



#separating the variables
train_text = train_text_org %>%
  separate(text, into = c("ID", "txt"), sep = "\\|\\|")
#str(train_text)

#id should be integer
train_text = train_text %>%
  mutate(ID = as.integer(ID))

#test text
test_text_org = tibble(text = read_lines('test_text', skip = 1))
test_text = test_text_org %>%
  separate(text, into = c("ID", "txt"), sep = "\\|\\|")

#ID should be integer
test_text = test_text %>%
  mutate(ID = as.integer(ID))

rm(train_text_org,test_text_org)


#The variants data tables
#train variants
train_data = train_data %>%
  mutate(Gene = factor(Gene),
         Variation = factor(Variation),
         Class = factor(Class))
#train_data$Class=as.factor(train_data$Class)
#text variants
test_data = test_data %>%
  mutate(Gene = factor(Gene),
         Variation = factor(Variation))

summary(train_data, maxsum = 9)
str(train)
str(test)

#checking for na's 
sum(is.na(train_data))
sum(is.na(test_data))

#total number of classes
levels(train$Class)

#summarizing the train/test data by Gene
#train variants
sum_train_g=train_data %>%
  group_by(Gene) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

sum_test_g=test_data %>%
  group_by(Gene) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

#summarizing the train/test data by variation
sum_train_v=train_data %>%
  group_by(Variation) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

sum_test_v=test_data %>%
  group_by(Variation) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

#Individual feature visualisations
#This is the frequency distribution of the most frequent Gene values:
#For train
top_gene = train_data %>%
  group_by(Gene) %>%
  summarise(count = n()) %>%
  filter(count > 50) %>% arrange(desc(count))
  


top_gene %>%
  ggplot(aes(reorder(Gene, -count, FUN = min), count)) +
  geom_point(size = 4,colour="orange") +
  labs(x = "Gene", y = "Frequency") + 
  coord_flip()


#For test
top_gene_test = test_data %>%
  group_by(Gene) %>%
  summarise(count = n()) %>%
  filter(count > 30) %>% arrange(desc(count))

top_gene_test %>%
  ggplot(aes(reorder(Gene, -count, FUN = min), count)) +
  geom_point(size = 4,colour="blue") +
  labs(x = "Gene", y = "Frequency") +
  coord_flip()

data_set_train = train_data %>% mutate(set = factor("train")) %>% select(-Class, -ID)
data_set_test = test_data %>% mutate(set = factor("test")) %>% select(-ID)

data_set_train = full_join(data_set_train, data_set_test)
rm(data_set_test)

#These are the most frequent Variations in the train (blue) vs test (red) data
data_set_train %>%
  group_by(Variation, set) %>%
  summarise(count = n()) %>%
  filter(count > 3) %>%
  ggplot(aes(reorder(Variation, -count, FUN = median), count, colour = set)) +
  geom_point(size = 4) + 
  coord_cartesian(ylim = c(0, 100)) + 
  labs(x = "Variation", y = "Frequency") 

# Class target is distribution in the train data:
train_data %>%
  ggplot(aes(Class)) +
  geom_bar(fill="pink")

#Feature interactions
#Gene vs Class
#logarithmic frequency scale to normalize the frequency
train_data %>%
  filter(Gene %in% str_c(top_gene$Gene)) %>%
  ggplot(aes(Gene)) +
  geom_bar(fill="red") +
  scale_y_log10() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=7,color = "blue"))+
  facet_wrap(~ Class)


# Classes sorted by Genes 
sorted_class_gen=train_data %>%
  filter(Gene %in% str_c(top_gene$Gene)) %>%
  ggplot(aes(Class)) +
  geom_bar(fill="orange") +
  scale_y_log10() +
  facet_wrap(~ Gene)

#Gene vs Variation
#for train
data_set_train = train_data %>%
  filter(Gene %in% str_c(top_gene$Gene)) %>%
  group_by(Gene, Variation) %>%
  summarise(count = n())

y_labels = str_sub(data_set_train$Variation, start = 1, end = 5)

data_set_train %>%
  ggplot(aes(reorder(Gene, count, FUN = median), reorder(Variation, count, FUN = median),colour = "orange")) +
  geom_count() +
  labs(x = "Gene", y = "Variation") +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=7),
        axis.ticks = element_blank(), axis.text.y = element_blank(),
        legend.position = "none")

#Test data
data_set_test = test_data %>%
  filter(Gene %in% str_c(top_gene$Gene)) %>%
  group_by(Gene, Variation) %>%
  summarise(count = n())

y_labels = str_sub(foo$Variation, start = 1, end = 5)

data_set_test %>%
  ggplot(aes(reorder(Gene, count, FUN = median), reorder(Variation, count, FUN = median),colour = "orange")) +
  geom_count() +
  labs(x = "Gene", y = "Variation") +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=7),
        axis.ticks = element_blank(), axis.text.y = element_blank(),
        legend.position = "none")


#The text files
text_samp = str_sub(train_text$txt[1], start = 1, end = 500)

#Feature Engineering
#Text length 
#for train
train_text = train_text %>%
  mutate(text_length = str_length(txt),
         set = "train")

#for test
test_text = test_text %>%
  mutate(text_length = str_length(txt),
         set = "test")

merge_text = full_join(train_text,test_text)

#overall distribution of the text entry lengths in train vs test

merge_text %>%
ggplot(aes(text_length, fill = set)) +
geom_histogram(bins = 50) + scale_fill_manual("SET", values = c("train" = "blue", "test" ="orange")) +  
labs(x = "Length of text")


#Facet wrap comparison of distribution changes for the different target Classes:
data_set_train = train_text %>%
  select(ID, text_length)
data_set_test <- train_data %>%
  select(ID, Class)


full_join(data_set_train, data_set_test, by = "ID") %>%
  ggplot(aes(text_length)) +
  geom_density(fill = "orange", bw = 5e3) +  
  labs(x = "Length of text") +
  facet_wrap(~ Class)



#empirical cumulative density functions
full_join(data_set_train, data_set_test, by = "ID") %>%
  ggplot(aes(text_length)) +
  stat_ecdf(geom = "step") +
  stat_ecdf(aes(text_length, color = Class), geom = "step") +
  labs(x = "Length of text",y="emperical_cummulative_density")

# median lengths for each class
data_set_train = train_text %>%
  select(ID, text_length)
bar = train_data %>%
  select(ID, Class)

avg_med = full_join(data_set_train, data_set_test, by = "ID") %>%
  group_by(Class) %>%
  summarise(med_text_length = median(text_length))

#bar plot for the result
ggplot(data=avg_med,aes(x=Class,y=med_text_length)) +
  geom_bar(stat="identity",fill="dark green")

#missing values
missing_values = merge_text %>%
  filter(text_length < 10)

#Keyword frequency - pedestrian approach
train_text = train_text %>%
  mutate( benign_freq = str_count(txt, "benign"),
          pathogenic_freq = str_count(txt, "pathogenic"))
        
  
#The frequency distributions of the word "pathogenic" for our 9 classes
#used logarithmic scale
  
data_set_train = train_text %>%
  select(ID, benign_freq, pathogenic_freq)
bar = train_data %>%
  select(ID, Class)

full_join(data_set_train, data_set_test, by = "ID") %>%
  ggplot(aes(pathogenic_freq,colour="red")) +
  geom_bar(fill="red") +
  scale_y_log10() +
  facet_wrap(~ Class)

# plot the ratio of the mean occurence per class of the word 
#"pathogenic" over the mean occurence of the word "benign":
data_set_train = train_text %>%
  select(ID, benign_freq, pathogenic_freq)
data_set_test = train_data %>%
  select(ID, Class)

full_join(data_set_train, data_set_test, by = "ID") %>%
  group_by(Class) %>%
  summarise(mean_benign = mean(benign_freq),
            mean_pathogenic = mean(pathogenic_freq),
            path_ben = mean(pathogenic_freq)/mean(benign_freq)) %>%
  ggplot(aes(reorder(Class, -path_ben, FUN = max), path_ben)) +
  geom_point(colour = "blue", size = 3) +
  labs(x = "Class", y = "pathogenic / benign")

#text analysis
#use the properties of tidytext
tdy = train_text %>% select(ID, txt) %>% unnest_tokens(word, txt)
head(tdy)

#defining list of own stopwords
data("stop_words")
my_stopwords = data_frame(word = c(as.character(1:100),
                                    "perform","discuss", "et", "al", "table","demonstrate",
                                    "data","method", "analysis", "analyze", "study",
                                     "fig", "figure","result", "conclusion", "author",
                                    "find", "found", "show", 
                                     "evaluate"))
tdy = tdy %>%
  anti_join(stop_words, by = "word") %>%
  anti_join(my_stopwords, by = "word") %>%
  filter(str_detect(word, "[a-z]"))

#plot for most occurance words
tdy %>%
  count(word) %>%
  filter(n > 6e4) %>%
  mutate(word = reorder(word, n)) %>%
  ggplot(aes(word,n)) +
  geom_col(fill = "orange") +
  xlab("words") + ylab("frequency") +
  coord_flip()

#after root stemming
# "SnowballC" package and its "wordStem" tool
tdy = tdy %>%
  mutate(word = wordStem(word))

tdy %>%
  count(word) %>%
  filter(n > 6e4) %>%
  mutate(word = reorder(word, n)) %>%
  ggplot(aes(word, n)) +
  geom_col(fill = "orange") +
  xlab("words") + ylab("frequency") +
  coord_flip()

#word cloud
#wordCloud library
pal2 = brewer.pal(8,"Dark2")
tdy %>% 
  count(word) %>%
  with(wordcloud(word, n, max.words = 100,random.order=F, colors=pal2))

#frequencies of words dependenting on some class
data_set_train = train_data %>%
  select(ID, Class)

tdy_class = full_join(tdy, data_set_train, by = "ID")

frequency = tdy_class %>%
  count(Class, word) %>%
  group_by(Class) %>%
  mutate(freq = n / sum(n)) %>% 
  select(-n) %>% 
  spread(Class, freq) %>% 
  gather(Class, freq, `1`:`2`)


#Class7 (the most frequent one) with Classes 1 and 2
#To get an idea taking only 500 occurences per Class to keep an overview
ggplot(frequency, aes(x = freq, y = `7`, color = abs(`7` - freq))) +
  geom_abline(color = "gray40", lty = 2) +
  geom_jitter(alpha = 0.1, size = 2.5, width = 0.1, height = 0.1) +
  geom_text(aes(label = word), check_overlap = TRUE, vjust = 1.5) +
  scale_x_log10(labels = percent_format()) +
  scale_y_log10(labels = percent_format()) +
  facet_wrap(~Class, ncol = 2) +
  theme(legend.position="none") 

#the correlation coefficients for each frequency for all the occurances
frequency = tdy_class %>%
  count(Class, word) %>%
  group_by(Class) %>%
  mutate(freq = n / sum(n)) %>% 
  select(-n) %>% 
  spread(Class, freq)

frequency %>%
  select(-word) %>%
  cor(use="complete.obs", method="spearman") %>%
  corrplot(type="lower", method="number", diag=FALSE)

#class 9 with class 3 and 5
data_set_train = train_data %>%
  select(ID, Class)

tyd_class = full_join(tdy, data_set_train, by = "ID")

frequency = tyd_class %>%
  count(Class, word) %>%
  filter(n > 20) %>%
  group_by(Class) %>%
  mutate(freq = n / sum(n)) %>% 
  select(-n) %>% 
  spread(Class, freq) %>% 
  gather(Class, freq, `3`,`5`)

#abline plot to compare class 9 with 3 & 5
ggplot(frequency, aes(x = freq, y = `9`, color = abs(`9` - freq))) +
geom_abline(color = "gray40", lty = 2) +
geom_jitter(alpha = 0.1, size = 2.5, width = 0.1, height = 0.1) +
  geom_text(aes(label = word), check_overlap = TRUE, vjust = 1.5) +
  scale_x_log10(labels = percent_format()) +
  scale_y_log10(labels = percent_format()) +
  facet_wrap(~Class, ncol = 2) +
  theme(legend.position="none") +
  labs(y = "Class 9", x = NULL)

#bind_tf_idf use to extract the metrics from a tidy data set that 
#contains words class and their counts/Class:

frequency = tdy_class %>%
  count(Class, word)

tf_idf = frequency %>%
  bind_tf_idf(word, Class, n)

#visualise  most characteristics words and their classes
tf_idf %>%
  arrange(desc(tf_idf)) %>%
  mutate(word = factor(word, levels = rev(unique(word)))) %>%
  top_n(20, tf_idf) %>%
  ggplot(aes(word, tf_idf, fill = Class)) + 
  geom_col() +scale_fill_manual("class", values = c("3" = "blue", "8" ="orange","9"="red")) +
  labs(x = NULL, y = "tf-idf") +
  coord_flip()

#overview of the most characteristic terms in each individual Class:
tf_idf %>%
  arrange(desc(tf_idf)) %>%
  mutate(word = factor(word, levels = rev(unique(word)))) %>%
  group_by(Class) %>%
  top_n(10, tf_idf) %>%
  ungroup() %>%  
  ggplot(aes(word, tf_idf, fill = Class)) +
  geom_col() + 
  labs(x = NULL, y = "tf-idf") +
  theme(legend.position = "none") +
  scale_fill_manual("class", values = c("1" = "orange", "2" ="sky blue","3"="red",
                                        "4" = "brown", "5" ="blue","6"="violet",
                                        "7" = "magenta", "8" ="dark green","9"="yellow")) +
  facet_wrap(~ Class, ncol = 3, scales = "free") +
  coord_flip()

##Word pair frequencies(n-grams concept)
tdy1 = train_text %>% select(ID, txt) %>% unnest_tokens(bigram, txt, token = "ngrams", n = 2)
head(tdy1)
tdy1 = tdy1[1:1000000,]

#to filter out the stopwords separate the bigrams first, 
#and then later unite them back together after the filtering. 
bi_sep = (tdy1 %>%
  separate(bigram, c("word1", "word2"), sep = " "))

bi_filt = bi_sep %>%
  filter(!word1 %in% stop_words$word) %>%
  filter(!word2 %in% stop_words$word) %>%
  filter(!word1 %in% my_stopwords$word) %>%
  filter(!word2 %in% my_stopwords$word)

#counting bi-gram combinations
bigram_counts = bi_filt %>%
  count(word1, word2, sort = TRUE)

#unite the words into bigram 
tdy1 = bi_filt %>%
  unite(bigram, word1, word2, sep = " ")

#Estimate tf-idf:
data_set_train = train_data %>%
  select(ID, Class)

tdy1_class = full_join(tdy1, foo, by = "ID")

#calculating tf-idf values
tdy1_tf_idf = tdy1_class %>%
  count(Class, bigram) %>%
  bind_tf_idf(bigram, Class, n) %>%
  arrange(desc(tf_idf))

#plot the bigrams per Class with the best tf-idf values:
tdy1_tf_idf %>%
  arrange(desc(tf_idf)) %>%
  mutate(bigram = factor(bigram, levels = rev(unique(bigram)))) %>%
  group_by(Class) %>%
  top_n(10, tf_idf) %>%
  ungroup() %>%  
  ggplot(aes(bigram, tf_idf, fill = Class)) +
  geom_col() +
  labs(x = NULL, y = "tf-idf") +
  theme(legend.position = "none") +
  facet_wrap(~ Class, ncol = 3, scales = "free") +
  coord_flip()
