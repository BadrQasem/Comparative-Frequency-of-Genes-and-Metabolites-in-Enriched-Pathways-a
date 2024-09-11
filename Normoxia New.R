attach(Metabolites_Normoxia)
head(Metabolites_Normoxia)
attach(Genes_Normoxia)
head(Genes_Normoxia)


# Check column names first to make sure we are using the correct names
names(Genes_Normoxia)
names(Metabolites_Normoxia)

# Now adjust renaming based on actual column names
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Adjust renaming based on correct column names
genes_data <- Genes_Normoxia %>%
  rename(Entity = Genes, Pathways = Pathways, Frequency = Frequency, Significance = Significance) %>%
  mutate(Type = "Gene")

# Use the actual name of the 'Metabolite' column here
metabolites_data <- Metabolites_Normoxia %>%
  rename(Entity = Metabolites, Pathways = Pathways, Frequency = Frequency, Significance = Significance) %>%
  mutate(Type = "Metabolite")


# Convert the Frequency column in both datasets to numeric
genes_data$Frequency <- as.numeric(genes_data$Frequency)
metabolites_data$Frequency <- as.numeric(metabolites_data$Frequency)

# Now combine the datasets
combined_data <- bind_rows(genes_data, metabolites_data)


# Combine the datasets
combined_data <- bind_rows(genes_data, metabolites_data)


aggregated_data <- combined_data %>%
  group_by(Pathways, Entity, Type) %>%
  summarise(Total_Frequency = sum(Frequency), .groups = 'drop')





# Separate heatmaps for Genes and Metabolites
library(gridExtra)

# Heatmap for Genes
heatmap_genes <- ggplot(aggregated_data %>% filter(Type == "Gene"), aes(x = Pathways, y = Entity, fill = Total_Frequency)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = c("lightblue", "blue"),
    values = scales::rescale(c(1, 2, 3, 4)),
    breaks = c(1, 2, 3, 4),
    labels = c("1", "2", "3", "4"),
    name = "Frequency (Genes)"
  ) +
  labs(title = "Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =  7),
        axis.text.y = element_text(size = 5),
        legend.position = "bottom")

# Heatmap for Metabolites
heatmap_metabolites <- ggplot(aggregated_data %>% filter(Type == "Metabolite"), aes(x = Pathways, y = Entity, fill = Total_Frequency)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = c("lightgreen", "darkgreen"),
    values = scales::rescale(c(1, 2, 3, 4, 5, 6, 7)),
    breaks = c(1, 2, 3, 4, 5, 6, 7),
    labels = c("1", "2", "3", "4", "5", "6", "7"),
    name = "Frequency (Metabolites)"
  ) +
  labs(title = "Metabolites") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom")

# Arrange the two plots side by side
grid.arrange(heatmap_genes, heatmap_metabolites, ncol = 2)












library(igraph)
library(ggraph)
library(dplyr)
library(ggplot2)

# Create the edge list (links between metabolites and pathways)
edges <- Metabolites_Normoxia %>%
  select(Metabolites, Pathways) %>%
  distinct()

# Calculate metabolite frequency
metabolites_freq <- edges %>%
  group_by(Metabolites) %>%
  summarise(frequency = n())  # Count how many pathways each metabolite is in

# Create a data frame with unique nodes and their type (Metabolite or Pathway)
nodes_df <- data.frame(
  name = unique(c(edges$Metabolites, edges$Pathways)),
  type = ifelse(unique(c(edges$Metabolites, edges$Pathways)) %in% edges$Metabolites, "Metabolites", "Pathway")
)

# Merge the frequency data with the nodes
nodes_df <- left_join(nodes_df, metabolites_freq, by = c("name" = "Metabolites"))

# Replace NA frequencies for pathways with 0 (since pathways don't have frequencies)
nodes_df$frequency[is.na(nodes_df$frequency)] <- 0

# Ensure frequency is numeric
nodes_df$frequency <- as.numeric(nodes_df$frequency)

# Create an igraph object with the vertices (nodes) and edges (connections)
network_graph <- graph_from_data_frame(edges, vertices = nodes_df, directed = FALSE)

# Plot the network with distinct colors and shapes for metabolites and pathways
ggraph(network_graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.6), edge_colour = "gray", show.legend = FALSE) +
  
  # Use color for metabolites based on frequency, fixed blue color for pathways, and remove size variation
  geom_node_point(aes(color = ifelse(type == "Metabolites", frequency, NA), 
                      shape = type), 
                  fill = ifelse(nodes_df$type == "Pathway", "blue", NA),  
                  size = ifelse(nodes_df$type == "Metabolites", 6, 4)) +  # Increase size of shapes) +
  
  # Add node labels
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  
  # Apply a gradient color scale for metabolites and a fixed blue color for pathways
  scale_color_gradient(low = "#ffeda0", high = "#f03b20", name = "Frequency of Metabolites", na.value = "blue") +
  
  # Assign fixed shapes to metabolites and pathways
  scale_shape_manual(values = c("Metabolites" = 16, "Pathway" = 17)) +  # 16 is a circle, 17 is a triangle
  
  
  # Increase shape size in the legend
  guides(shape = guide_legend(override.aes = list(size = 5))) +

  
  theme_void() +
  theme(legend.position = "right") +
  labs(title = "Network of Metabolites with Enriched Pathways")







library(igraph)
library(ggraph)
library(dplyr)
library(ggplot2)

# Create the edge list (links between Genes and pathways)
edges <- Genes_Hypoxia %>%
  select(Genes, Pathways) %>%
  distinct()

# Calculate Gene frequency
Genes_freq <- edges %>%
  group_by(Genes) %>%
  summarise(frequency = n())  # Count how many pathways each Gene is in

# Create a data frame with unique nodes and their type (Gene or Pathway)
nodes_df <- data.frame(
  name = unique(c(edges$Genes, edges$Pathways)),
  type = ifelse(unique(c(edges$Genes, edges$Pathways)) %in% edges$Genes, "Genes", "Pathway")
)

# Merge the frequency data with the nodes
nodes_df <- left_join(nodes_df, Genes_freq, by = c("name" = "Genes"))

# Replace NA frequencies for pathways with 0 (since pathways don't have frequencies)
nodes_df$frequency[is.na(nodes_df$frequency)] <- 0

# Ensure frequency is numeric
nodes_df$frequency <- as.numeric(nodes_df$frequency)

# Create an igraph object with the vertices (nodes) and edges (connections)
network_graph <- graph_from_data_frame(edges, vertices = nodes_df, directed = FALSE)

# Plot the network with distinct colors and shapes for Genes and pathways
ggraph(network_graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.6), edge_colour = "gray", show.legend = FALSE) +
  
  # Use color for Genes based on frequency, fixed blue color for pathways, and remove size variation
  geom_node_point(aes(color = ifelse(type == "Genes", frequency, NA), 
                      shape = type), 
                  fill = ifelse(nodes_df$type == "Pathway", "blue", NA),  
                  size = ifelse(nodes_df$type == "Genes", 6, 4)) +  # Increase size of shapes) +
  
  # Add node labels
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  
  # Apply a gradient color scale for Genes and a fixed blue color for pathways
  scale_color_gradient(low = "#ffeda0", high = "#f03b20", name = "Frequency of Genes", na.value = "blue") +
  
  # Assign fixed shapes to Genes and pathways
  scale_shape_manual(values = c("Genes" = 16, "Pathway" = 17)) +  # 16 is a circle, 17 is a triangle
  
  
  # Increase shape size in the legend
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  
  
  theme_void() +
  theme(legend.position = "right") +
  labs(title = "Network of Genes with Enriched Pathways")









