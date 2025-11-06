from nicegui import ui,app
from model import linear_chain
from content import background_markdown
from ui_helpers import section_header, divider, param_row, build_info_dialog, header_with_info, make_echart_card
from pathlib import Path
from model.parameters import default_params, drivers_from_params  # ← your real params & seasonal drivers
from model.fest_params_run_SH_out import runSH
# from plotting import plot_panel


# Serve the local "static" folder at the URL path "/static"
app.add_static_files('/static', str(Path(__file__).parent / 'static'))

ui.page_title('Coral DEB model')

info_dialog, open_info = build_info_dialog(background_markdown())
header_with_info('Coral DEB model', open_info)

with ui.element('div').classes('flex flex-col md:flex-row gap-6 p-6 w-full max-w-[1800px] mx-auto'):
    with ui.card().classes('w-full md:w-[420px] shrink-0 shadow-md p-6'):
        section_header('Parameters').classes('mb-4')

        n, _ = param_row('Number of sub-stages n', min=1, max=100, value=10, step=1, fmt=lambda v: f'{int(v)}')
        divider()
        mu, _ = param_row('Mean delay μ', min=1.0, max=8.0, value=3.0, step=0.1, fmt=lambda v: f'{v:.2f}')
        divider()
        section_header('Properties of the delay').classes('mb-2')
        mean_label = ui.label().classes('mt-1')
        var_label  = ui.label().classes('mt-1')

    chart = make_echart_card('S/H ratio')

def update_chart(_=None):
    mu_val = float(mu.value)
    n_val = int(n.value)
#    tmax_val = 10.0

    mean_label.text = f'mean = {mu_val:.4g}'
    var_label.text  = f'variance = {mu_val ** 2 / n_val:.4g}'
    parameters = default_params()
    parameters['T_amp']= mu_val
    t, pdf = runSH(parameters)#linear_chain(mu_val, n_val, tmax=tmax_val)
    tmax_val = t[-1]
    data = [[t[i], pdf[i]] for i in range(len(t))]

    # Update single series (combined line + area)
    chart.options['series'][0]['data'] = data
    chart.options['xAxis']['max'] = tmax_val
    chart.update()


for w in (mu, n):
    w.on('change', update_chart)

update_chart()
ui.add_head_html('''
<link rel="icon" type="image/x-icon" href="/static/favicon.ico?v=10">
''')
ui.run(host='0.0.0.0', port=8080)
