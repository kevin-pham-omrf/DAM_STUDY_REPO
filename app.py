import pyarrow.parquet as pq
import pandas as pd
import plotly.graph_objects as go
import numpy as np
from scipy.stats import norm
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shinywidgets import output_widget, render_widget

genes = pd.read_csv("DETECTED_GENES.csv", header=None, index_col=False)

mode="light"

population = ["DAM", "Old Homeostatic", "Transition", "Young Homeostatic"]
treatment = ["AD", "EAE", "Aging"]

uni_mean = -0.032989667
uni_sd = 0.289933333
sim_mean = -0.1818
sim_sd = 0.5658
sim_x_values = np.linspace(sim_mean - 4 * sim_sd, sim_mean + 4 * sim_sd, 1000)
sim_y_values = norm.pdf(sim_x_values, loc=sim_mean, scale=sim_sd)
uni_x_values = np.linspace(uni_mean - 4 * uni_sd, uni_mean + 4 * uni_sd, 1000)
uni_y_values = norm.pdf(uni_x_values, loc=uni_mean, scale=uni_sd)

def sort_key(g):
    starts_with_digit = g[0].isdigit()
    ends_with_rik = g.lower().endswith('rik')
    return (ends_with_rik, starts_with_digit, g.lower())

sorted_genes = sorted(genes[0], key=sort_key)

color_map = {
    'DAM': "#74a5ce",
    'Old Homeostatic': "#afc75b",
    'Transition': "#b14646",
    'Young Homeostatic': "#4c3ec7",
    'AD': "#b33ba9",
    'EAE': "#4d26bb",
    'Aging': "#53410F",
    'Old Homeostatic AD': "#bb58da",
    'DAM AD': "#4ecbeb",
    'Old Homeostatic EAE': "#f088e7",
    'DAM EAE': "#6d93e6",
    'Old Homeostatic Aging': "#c65cdb",
    'DAM Aging': "#4224ee",
    'Transition Aging': "#eb25eb",
    'Young Homeostatic Aging': "#9694f3"
}

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_selectize(
            "gene", 
            "Gene", 
            sorted_genes,
            selected="Cx3cr1"
        ),
        ui.input_select(
            "filter", 
            "Grouping", 
            {
                1: "Treatment",
                2: "Population",
                3: "Treatment and Population"
            },
            selected=1
        ),
        ui.panel_conditional(
            "input.filter == 1 || input.filter == 3",
            ui.input_checkbox_group(
                "treatment",
                "Filter by Treatment",
                {
                    "AD": "AD",
                    "EAE": "EAE",
                    "Aging": "Aging"
                },
                selected=treatment
            ),
        ),
        ui.panel_conditional(
            "input.filter == 2 || input.filter == 3",
            ui.input_checkbox_group(
                "population",
                "Filter by Population",
                {
                    "Young Homeostatic": "Young Homeostatic",
                    "Old Homeostatic": "Old Homeostatic",
                    "Transition": "Transition", 
                    "DAM": "DAM"
                },
                selected=population
            ),
        ),
        ui.download_button("download_expr", "Download Expression Data"),
        ui.input_action_button("toggle_dark", "Toggle Dark Mode")
    ),
    ui.page_navbar(
        ui.nav_panel(
            "Gene Expression",
            ui.layout_columns(
                output_widget("expression_plot")
            )
        ),
        ui.nav_panel(
            "Similarity and Uniqueness",
            ui.layout_columns(
                output_widget("similarity_plot")
            ),
            ui.layout_columns(
                output_widget("uniqueness_plot")
            )
        ),
    )
)

def server(input: Inputs, output: Outputs, session: Session):
    @render.download(filename="gene_expression_data.csv")
    def download_expr():
        outData = filtered_expr()
        yield outData.to_csv(index=False)
    
    @reactive.effect
    @reactive.event(input.toggle_dark)
    def _():
        global mode
        if mode == "light":
            ui.update_dark_mode("dark")
            mode = "dark"
        else:
            ui.update_dark_mode("light")
            mode = "light"
    
    @reactive.Calc
    def filtered_expr() -> pd.DataFrame:
        data = pq.read_table('ALL_RPKM_LABELED_TRANSPOSED.parquet', columns=["Treatment", "Subtype", input.gene()]).to_pandas()
        data["Subtype"] = pd.Categorical(data["Subtype"], categories=population, ordered=True)
        data["Treatment"] = pd.Categorical(data["Treatment"], categories=treatment, ordered=True)
        sorted_data = data.sort_values(['Treatment', 'Subtype'])
        match input.filter():
            case "1":
                outputData = sorted_data.loc[sorted_data["Treatment"].isin(input.treatment()), ("Treatment", input.gene())]
                outputData["GROUP"] = outputData["Treatment"]
                return outputData
            case "2":
                outputData = sorted_data.loc[sorted_data["Subtype"].isin(input.population()), ("Subtype", input.gene())]
                outputData["GROUP"] = outputData["Subtype"]
                return outputData
            case "3":
                outputData = sorted_data.loc[sorted_data["Subtype"].isin(input.population()), ("Subtype", "Treatment", input.gene())]
                outputData = outputData.loc[outputData["Treatment"].isin(input.treatment()), ]
                outputData["GROUP"] = outputData["Subtype"].astype(str) + " " + outputData["Treatment"].astype(str)
                return outputData

    @reactive.Calc
    def filtered_sim() -> pd.DataFrame:
        data = pq.read_table('ALL_SIM_SCORES_TRANSPOSED.parquet', columns=["Gene Symbol", input.gene()]).to_pandas()
        outputData = data.loc[data["Gene Symbol"].isin(["ALL_X", "ALL_Y"]), ("Gene Symbol", input.gene())]
        outputData = outputData.T
        outputData.columns = outputData.iloc[0]
        outputData = outputData[1:].reset_index(drop=True)
        return outputData

    @reactive.Calc
    def filtered_uni() -> pd.DataFrame:
        data = pq.read_table('ALL_SIM_SCORES_TRANSPOSED.parquet', columns=["Gene Symbol", input.gene()]).to_pandas()
        outputData = data.loc[data["Gene Symbol"].isin(["AD_X", "EAE_X", "AGE_X", "AD_Y", "EAE_Y", "AGE_Y"]), ("Gene Symbol", input.gene())]
        outputData = outputData.T
        outputData.columns = outputData.iloc[0]
        outputData = outputData[1:].reset_index(drop=True)
        return outputData

    @render_widget
    def expression_plot():
        if input.gene() not in sorted_genes:
            return
        fig = go.Figure()
        data = filtered_expr()

        if mode == "light":
            bgcolor = "#e4e4e4"
            fontcolor = "#000000"
            papercolor = "#ffffff"
        else:
            bgcolor = "#B9B9B9"
            fontcolor = "#ffffff"
            papercolor = "#252525"
        
        for group in data['GROUP'].unique():
            fig.add_trace(go.Box(y = data[data['GROUP'] == group][input.gene()],
                boxpoints = 'all', jitter = 0.5, marker_line_width=1, line = dict(width=2),
                pointpos = 0, name = group, marker_color=color_map[group]))
        fig.update_layout(
            title=input.gene(),
            title_font = dict(
                size = 24,
                textcase = "upper",
                weight = "bold",
                color = fontcolor
            ),
            font = dict(
                color = fontcolor
            ),
            yaxis_title = "RPKM",
            xaxis_title = "",
            showlegend = False,
            paper_bgcolor = papercolor,
            plot_bgcolor = bgcolor,
            margin=dict(t=100)
        )
        return fig
    
    @render_widget
    def similarity_plot():
        if input.gene() not in sorted_genes:
            return
        fig = go.Figure()
        data = filtered_sim()

        if mode == "light":
            bgcolor = "#e4e4e4"
            fontcolor = "#000000"
            papercolor = "#ffffff"
        else:
            bgcolor = "#B9B9B9"
            fontcolor = "#ffffff"
            papercolor = "#252525"
        
        fig.add_trace(go.Scatter(
            x = sim_x_values,
            y = sim_y_values,
            mode='lines'
        ))
        fig.add_trace(go.Scatter(
            x = data["ALL_X"],
            y = data["ALL_Y"],
            mode='markers',
            marker=dict(
                size=10,
                color='red',
                symbol='circle'
            )
        ))
        fig.update_layout(
            title="Similarity Score",
            title_font = dict(
                size = 24,
                textcase = "upper",
                weight = "bold",
                color = fontcolor
            ),
            font = dict(
                color = fontcolor
            ),
            yaxis_title = "",
            xaxis_title = "",
            showlegend = False,
            paper_bgcolor = papercolor,
            plot_bgcolor = bgcolor,
            margin=dict(t=100)
        )
        return fig
    
    @render_widget
    def uniqueness_plot():
        if input.gene() not in sorted_genes:
            return
        fig = go.Figure()
        data = filtered_uni()

        if mode == "light":
            bgcolor = "#e4e4e4"
            fontcolor = "#000000"
            papercolor = "#ffffff"
        else:
            bgcolor = "#B9B9B9"
            fontcolor = "#ffffff"
            papercolor = "#252525"
        
        fig.add_trace(go.Scatter(
            x = uni_x_values,
            y = uni_y_values,
            mode='lines'
        ))
        fig.add_trace(go.Scatter(
            x = [data["AD_X"].iloc[0]],
            y = [norm.pdf(data["AD_X"].iloc[0], loc = uni_mean, scale = uni_sd)],
            mode='markers',
            marker=dict(
                size=10,
                color='red',
                symbol='circle'
            )
        ))
        fig.add_trace(go.Scatter(
            x = [data["EAE_X"].iloc[0]],
            y = [norm.pdf(data["EAE_X"].iloc[0], loc = uni_mean, scale = uni_sd)],
            mode='markers',
            marker=dict(
                size=10,
                color='green',
                symbol='circle'
            )
        ))
        fig.add_trace(go.Scatter(
            x = [data["AGE_X"].iloc[0]],
            y = [norm.pdf(data["AGE_X"].iloc[0], loc = uni_mean, scale = uni_sd)],
            mode='markers',
            marker=dict(
                size=10,
                color='blue',
                symbol='circle'
            )
        ))
        fig.update_layout(
            title="Uniqueness Score",
            title_font = dict(
                size = 24,
                textcase = "upper",
                weight = "bold",
                color = fontcolor
            ),
            font = dict(
                color = fontcolor
            ),
            yaxis_title = "",
            xaxis_title = "",
            showlegend = False,
            paper_bgcolor = papercolor,
            plot_bgcolor = bgcolor,
            margin=dict(t=100)
        )
        return fig

app = App(app_ui, server)