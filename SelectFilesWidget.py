'''
Adapted from: https://gist.github.com/draperjames/90403fd013e332e4a14070aab3e3e7b0w
'''
import ipywidgets as widgets
from tkinter import Tk, filedialog
import traitlets

class SelectFilesButton(widgets.Button):
    """A file widget that leverages tkinter.filedialog."""

    def __init__(self):
        super(SelectFilesButton, self).__init__()
        # Add the selected_files trait
        self.add_traits(files=traitlets.traitlets.List())
        # Create the button.
        self.description = "Select file"
        self.icon = "square-o"
        self.style.button_color = "red"
        # Set on click behavior.
        self.on_click(self.select_files)
        self.layout.width='20%'

    @staticmethod
    def select_files(b):
        """Generate instance of tkinter.filedialog.

        Parameters
        ----------
        b : obj:
            An instance of ipywidgets.widgets.Button
        """
        # Create Tk root
        root = Tk()
        # Hide the main window
        root.withdraw()
        # Raise the root to the top of all windows.
        root.call('wm', 'attributes', '.', '-topmost', True)
        # List of selected fileswill be set to b.value
        ftypes = [('XRD data in xy format', '*.txt')]
        b.files = filedialog.askopenfilename(filetypes = ftypes, multiple=True)

        #b.description = "Files Selected"
        b.icon = "check-square-o"
        b.style.button_color = "lightgreen"
