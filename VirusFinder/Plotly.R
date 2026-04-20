##When importing dataset, remember to convert numerals to numric
library(plotly)
p <- plot_ly(data = Virus_SummaryCount_DotPlotBarGraph) %>% 
  add_trace(x = ~ICTV_Family_Name, 
        y = ~LTNP,
        name = 'Virus Species Count In LTNPs',
        type = 'scatter',
        mode = 'markers',
        opacity = 1,
        size = ~LTNP_SIZE,
        marker = list(
          color = c('#FFFFFF'),line = list(width = 0.5, color = '#000000')
        ))%>%
  add_trace(x = ~ICTV_Family_Name, 
            y = ~RP,
            name = 'Virus Species Count In RPs',
            type = 'scatter',
            mode = 'markers',
            opacity = 0.8,
            size = ~RP_SIZE,
            marker = list(
              color = c('#FFCF79'), line = list(width = 0.5, color = '#000000')
            ))%>%
  add_trace(x = ~ICTV_Family_Name, y = ~LTNP, name = 'Cumulative Count In LTNPs', type = 'bar',
            opacity = 1, marker = list(color = '#2C6700'))%>%
  add_trace(x = ~ICTV_Family_Name, y = ~RP, name = 'Cumulative Count In RPs', type = 'bar',
          opacity = 1, marker = list(color = '#92CD00'))%>%
  layout(title = '',
         xaxis = list(
           title = "Virus Family",
           titlefont = list(
             size = 18,
             color = 'rgb(0, 0, 0)'),
           tickfont = list(
             size = 14,
             color = 'rgb(0, 0, 0)')),
         yaxis = list(
           title = 'Number of Samples',
           titlefont = list(
             size = 18,
             color = 'rgb(0, 0, 0)'),
           tickfont = list(
             size = 14,
             color = 'rgb(0, 0, 0)')),
         legend = list(x = 0.65, y = 1, font = list(family = "sans-serif", size = 14, color = 'rgb(0, 0, 0)'), bgcolor = 'rgba(255, 255, 255, 0)', bordercolor = 'rgba(255, 255, 255, 0)'),
         barmode = 'group', bargap = 0, bargroupgap = 0)

p
###export as pdf
#install orca: https://github.com/plotly/orca#installation
if (!require("processx")) install.packages("processx")

orca(p, "disease_progression.svg")


p2 <- plot_ly(data = Virus_SummaryCount_DotPlotBarGraph) %>% 
  add_trace(x = ~ICTV_Family_Name, 
            y = ~BWR,
            name = 'Virus Species Count In BWRs',
            type = 'scatter',
            mode = 'markers',
            opacity = 1,
            size = ~BWR_SIZE,
            marker = list(color = c('#FFDE00'), line = list(width = 0.5, color = '#000000')
                          ))%>%
  add_trace(x = ~ICTV_Family_Name, 
            y = ~UGR,
            name = 'Virus Species Count In UGRs',
            type = 'scatter',
            mode = 'markers',
            opacity = 0.8,
            size = ~UGR_SIZE,
            marker = list(color = c('#097054'), line = list(width = 0.5, color = '#000000')
            ))%>%
  add_trace(x = ~ICTV_Family_Name, y = ~BWR, name = 'Cumulative Count In BWRs', type = 'bar',
            opacity = 1, marker = list(color = '#FF9900'))%>%
  add_trace(x = ~ICTV_Family_Name, y = ~UGR, name = 'Cumulative Count In UGRs', type = 'bar',
            opacity = 1, marker = list(color = '#6599FF'))%>%
  layout(title = '',
         xaxis = list(
           title = "Virus Family",
           titlefont = list(
             size = 18,
             color = 'rgb(0, 0, 0)'),
           tickfont = list(
             size = 14,
             color = 'rgb(0, 0, 0)')),
         yaxis = list(
           title = 'Number of Samples',
           titlefont = list(
             size = 18,
             color = 'rgb(0, 0, 0)'),
           tickfont = list(
             size = 14,
             color = 'rgb(0, 0, 0)')),
         legend = list(x = 0.65, y = 1, font = list(family = "sans-serif", size = 14, color = 'rgb(0, 0, 0)'), bgcolor = 'rgba(255, 255, 255, 0)', bordercolor = 'rgba(255, 255, 255, 0)'),
         barmode = 'group', bargap = 0, bargroupgap = 0)
p2
if (!require("processx")) install.packages("processx")

orca(p2, "country.svg")