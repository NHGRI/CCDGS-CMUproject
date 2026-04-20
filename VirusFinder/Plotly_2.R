###THIS IS FOR TESTING!!
Virus_SummaryCount_DotPlotBarGraph <- read_excel("C:/Users/savan/OneDrive/Cafgen/Analysis/Viral_Reads/Virus_SummaryCount_DotPlotBarGraph.xlsx", 
                                                   +     sheet = "Virus_SummaryCount_PROGRESSION")
View(Virus_SummaryCount_DotPlotBarGraph)
fig <- plot_ly(data = Virus_SummaryCount_DotPlotBarGraph) %>%
  add_trace(y = ~Count,
            type ='box',
            boxpoints = "all",
            jitter = 0.5,
            pointpos = -1.8,
            color = ~ICTV_Family_Name)%>%
  add_trace(data = Virus_SummaryCount_DotPlotBarGraph, x = ~ICTV_Family_Name, y = ~LTNP_RP, type = 'bar')


fig
 

