# Load necessary libraries
library(ggplot2)

# Sample data
data <- data.frame(
  Family = c("Family 1", "Family 1", "Family 1", "Family 2", "Family 2", 
             "Family 3", "Family 3", "Family 4", "Family 5"),
  Porphyria_Type = c("AIP", "AIP", "AIP", "VP", "VP", 
                     "EPP", "EPP", "fPCT", "CEP"),
  Patient_ID = c("P001", "P002", "P003", "P004", "P005", 
                 "P006", "P007", "P008", "P009")
)

# Summarize data to get family size
family_data <- data %>%
  group_by(Family, Porphyria_Type) %>%
  summarise(Family_Size = n()) %>%
  ungroup()

# Create the bubble chart
ggplot(family_data, aes(x = Family, y = Family_Size, size = Family_Size, fill = Porphyria_Type)) +
  geom_point(shape = 21, color = "black") +
  scale_size_area(max_size = 15) +  # Adjust the size of bubbles
  labs(title = "Porphyria Families: Size and Type", 
       x = "Family", 
       y = "Family Size") +
  theme_minimal() +
  theme(legend.title = element_blank())  # Remove the legend title for the fill color
