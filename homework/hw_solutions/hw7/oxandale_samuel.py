import sys
import platform
import numpy as np

from PyQt5.QtWidgets import (QMainWindow, QApplication, QDialog, QLineEdit, 
                             QVBoxLayout, QAction, QMessageBox,QFileDialog,
                             QSizePolicy, QComboBox)
from PyQt5.QtCore import QT_VERSION_STR, PYQT_VERSION_STR

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

class MainWindow(QMainWindow) :
    
    def __init__(self, parent=None) :
        QMainWindow.__init__(self, parent)

        self.setWindowTitle('Function Plotter')

        ########################################################################
        # ADD MENU ITEMS
        ########################################################################
        
        self.menuBar().setNativeMenuBar(False)
        
        # Create the File menu
        self.menuFile = self.menuBar().addMenu("&File")
        self.actionSaveAs = QAction("&Save As", self)
        self.actionSaveAs.triggered.connect(self.saveas)
        self.actionQuit = QAction("&Quit", self)
        self.actionQuit.triggered.connect(self.close)
        self.menuFile.addActions([self.actionSaveAs, self.actionQuit])
        
        # Create the Help menu
        self.menuHelp = self.menuBar().addMenu("&Help")
        self.actionAbout = QAction("&About",self)
        self.actionAbout.triggered.connect(self.about)
        self.menuHelp.addActions([self.actionAbout])
        
        ########################################################################
        # CREATE CENTRAL WIDGET
        ########################################################################

        self.widget = QDialog()
        
        self.plot = MatplotlibCanvas()

        self.cbx_equations = EquationCombo()
        self.cbx_equations.activated.connect(self.cbx_equations_updated)
        
        self.xvalues = QLineEdit("np.linspace(0,10.0)")
        self.xvalues.returnPressed.connect(self.update_plot)
        
        self.enter_equation = QLineEdit("enter f(x)")
        self.enter_equation.returnPressed.connect(self.update_cbx)

        layout = QVBoxLayout()
        layout.addWidget(self.plot)
        layout.addWidget(self.cbx_equations)
        layout.addWidget(self.xvalues) 
        self.widget.setLayout(layout)        
        self.setCentralWidget(self.widget) 

    def saveas(self):
        """Save x and f(x) to a *.txt file"""
        filename, filetype = QFileDialog.getSaveFileName(self,
                            "QFileDialog.getSaveFileName()", "",
                            "Text Files (*.txt)")
        eqn = str(self.cbx_equations.currentText())
        x = eval(self.xvalues.text())
        y = eval(eqn)
        np.savetxt(filename, np.array([x,y]).T)
 
        
    def about(self):
        QMessageBox.about(self, 
            "About Fucnction Plotter",
            """<b>Function Plotter</b>
               <p>Copyright &copy; 2017 Samuel Oxandale, All Rights Reserved.
               <p>Python %s -- Qt %s -- PyQt %s on %s""" %
            (platform.python_version(),
             QT_VERSION_STR, PYQT_VERSION_STR, platform.system()))
        
    def cbx_equations_updated(self) :
        if str(self.cbx_equations.currentText()) == "enter custom equation" :
            self.custom_equation()
        else :
            self.update_plot()
            
    def update_cbx(self) :
        custom_eq = str(self.enter_equation.text())
        self.cbx_equations.add_select_item(custom_eq)
        self.update_plot()
        self.enter_equation.hide()

    def custom_equation(self) :
        self.enter_equation.show()

    def update_plot(self) : 
        eqn = str(self.cbx_equations.currentText())
        x = eval(self.xvalues.text())
        y = eval(eqn)
        self.plot.redraw(x,y)
        
class MatplotlibCanvas(FigureCanvas) :
    """ This is borrowed heavily from the matplotlib documentation;
        specifically, see:
        http://matplotlib.org/examples/user_interfaces/embedding_in_qt5.html
    """
    def __init__(self):
        
        # Initialize the figure and axes
        self.fig = Figure()
        self.axes = self.fig.add_subplot(111)
        
        # Give it some default plot (if desired).  
        x = np.arange(0.0, 10.0, 0.01)
        y = np.sin(x)
        self.axes.plot(x, y)
        self.axes.set_xlabel('x')
        self.axes.set_ylabel('f(x)')  
        self.axes.set_title('f(x) vs. x')
        
        # Now do the initialization of the super class
        FigureCanvas.__init__(self, self.fig)
        #self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
         
        
    def redraw(self, x, y) :
        """ Redraw the figure with new x and y values.
        """
        # clear the old image (axes.hold is deprecated)
        self.axes.clear()
        self.axes.plot(x, y)
        self.axes.set_xlabel('x')
        self.axes.set_ylabel('f(x)')  
        self.axes.set_title('f(x) vs. x')
        self.axes.plot(x, y)
        self.draw()    
        
    def set_title(self, s):
        self.axes.set_title(s)
        self.draw()  
        
class EquationCombo(QComboBox) :
    
    def __init__(self) :
        QComboBox.__init__(self)
        # Set items in ComboBox
        equation_list = ["np.sin(x)", "np.cos(x)", "np.exp(x)"]
        for eq in equation_list :
            self.addItem(eq)
        self.addItem("enter custom equation", 99)
        
    def add_select_item(self, item) :
        self.addItem(item)
        index = self.findText(item)
        if index >= 0:
            self.setCurrentIndex(index)
        
    def new_selection(self) :
        pass
    
        
        
if __name__ == "__main__" :        
    app = QApplication(sys.argv)
    widget = MainWindow()
    widget.show()
    app.exec_()
