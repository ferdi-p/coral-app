from nicegui import ui
from content import background_markdown
from pathlib import Path
from nicegui import app

# serve your static folder so images like /static/diagram.svg work
app.add_static_files('/static', str(Path(__file__).parent / 'static'))

# show the background markdown with LaTeX enabled
ui.markdown(background_markdown(), extras=['latex']).classes('prose mx-auto p-6 max-w-[900px]')

# optional: tweak styling
ui.add_head_html('''
<style>
  body { background-color: #fafafa; }
  .prose p:has(img) { margin-top: 0 !important; margin-bottom: 0 !important; }
</style>
''')

ui.run(title='Background Preview', host='0.0.0.0', port=8080, show=True)
