from nicegui import ui
from typing import Callable, Tuple, Union

def section_header(text: str):
    return ui.label(text).classes('text-lg font-medium')

def divider(margin: str = 'my-5', opacity: str = 'opacity-60'):
    return ui.separator().classes(f'{margin} {opacity}')

def param_row(
    title: str,
    *,
    min: Union[float, int],
    max: Union[float, int],
    value: Union[float, int],
    step: Union[float, int],
    fmt: Callable[[float], str],
) -> Tuple[ui.slider, ui.label]:
    ui.label(title).classes('font-medium mb-1')
    with ui.element('div').classes('flex items-center gap-4 w-full'):
        s = ui.slider(min=min, max=max, value=value, step=step).classes('flex-1')
        v = ui.label(fmt(float(value))).classes('w-16 text-right')
        v.bind_text_from(s, 'value', backward=lambda x: fmt(float(x)))
    return s, v

def build_info_dialog(markdown_text: str):
    with ui.dialog() as dialog, ui.card().classes('w-[min(95vw,900px)] max-h-[80vh] overflow-auto p-6'):
        ui.markdown(markdown_text, extras=['latex'])
        ui.button('Close', on_click=dialog.close).classes('mt-3')
    def open_dialog():
        dialog.open()
    return dialog, open_dialog

def header_with_info(title: str, on_info_click):
    # Title and icon button side by side on the left
    with ui.header().classes('bg-gray-200 shadow-sm'):
        with ui.row().classes('items-center gap-6 px-5 py-3'):
            ui.label(title).classes('text-xl font-semibold text-gray-900')
            ui.button(icon='info', on_click=on_info_click).props('round unelevated size=md color=primary')


def make_echart_card(title: str, height: str = '50vh', max_width: str = '800px'):
    # Card expands to fill the right side of the flex row
    with ui.card().classes('w-full md:w-[800px] shrink-0 shadow-md p-4'):
        ui.label(title).classes('text-lg font-medium mb-2 text-center')
        # Centered width anchor; chart will fill this
        with ui.element('div').style(f'width:100%; max-width:{max_width}; margin:0 auto;'):
            chart = ui.echart({
                'tooltip': {'show': False},
                'grid': {'left': 60, 'right': 30, 'top': 30, 'bottom': 40},
                'xAxis': {'type': 'value', 'name': 't', 'min': 0},
                'yAxis': {'type': 'value', 'name': 'density'},
                'series': [
                    {'type': 'line', 'showSymbol': False, 'areaStyle': {}, 'data': []},
                ],
            }).style(f'height:{height}; width:100%')  # fill the wrapper
            return chart
