from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QHBoxLayout,
                             QLabel, QSizePolicy, QSlider, QSpacerItem,
                             QVBoxLayout, QWidget, QGridLayout)

from PyQt5 import QtWidgets, QtCore, QtGui
import pyqtgraph as pg
import numpy as np
import sys
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
# Import pandas package 
import pandas as pd


devtools = importr('devtools')
devtools.install_github('https://github.com/mankeldy/Methanogen_Package')
stats = importr('stats')
grdevices = importr('grDevices')
base = importr('base')
datasets = importr('datasets')
CHNOSZ = importr('CHNOSZ')
microbialkitchen=importr('microbialkitchen',on_conflict="warn")
methanogen = importr('Methanogen')

colors =  {
            'lightest':"#eeeeee",
            'lighter':"#e5e5e5",
            'light':"#effffb",
            'himid':"#50d890",
            'midmid':"#1089ff",
            'lomid':"#4f98ca",
            'dark' :"#272727",
            'darker' :"#23374d",
}

pg.setConfigOptions(antialias=True)
pg.setConfigOption('background', colors['dark'])
pg.setConfigOption('foreground', colors['light'])

QLabel_style = f"""
QLabel{{
    color: {colors['light']};
}}
"""


def fac(CH4,H2,DIC,pH,temperature,volumesoln,volumehead):
    output = methanogen.methanogenesis(CH4_initial = CH4*1e-6,
                                       H2_initial = 1e-3,
                                       DIC_initial = DIC*1e-3,
                                       pH_initial = pH,
                                       standard_gibbs = -191359.46584,
                                       temperature = temperature + 273.15,
                                       VolumeSolution = volumesoln*1e-3,
                                       VolumeHeadspace = volumehead*1e-3,
                                       delta_DIC = 0.0001, 
                                       biomass_yield = 2.4,
                                       carbon_fraction=0.44)
    
    with localconverter(robjects.default_converter + pandas2ri.converter):
        pd_from_r_df = robjects.conversion.rpy2py(output)
    return pd_from_r_df
print("before")
names = list(fac(1,1,1,1,1,1,1).columns)
print("after")
def get_vals(function,y,**kwargs):
    
    run = function(**kwargs)
    x_vals = run['DIC.consumed']
    y_vals = run[y]
    return x_vals,y_vals

def scatter(plt,x_vals,y_vals):
    plt.clear()
    line = plt.plot(x=x_vals, y=y_vals)
    line.setPen(color=colors['lomid'], width=3.0)

    scat = pg.ScatterPlotItem(x=x_vals, y=y_vals)
    scat.setPen(color=(1,1,1,0), width=0.0)
    plt.addItem(scat)
    #plt.setYRange(0,14)
    #plt.setXRange(0,1e-2)

def linear_interp(x, in_min, in_max, out_min, out_max, lim=True, dec=None):
    y = (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min

    if dec is not None:
        y = round(y, dec)

    if not lim:
        return y

    else:
        y = max(out_min, y)
        y = min(out_max, y)

        return y

class Controls(QWidget):
    def __init__(self, variable='', parent=None):
        super(Controls, self).__init__(parent=parent)
        
        self.controlpanelLayout = QGridLayout(self)
        self.environmentLayout = QVBoxLayout(self)
        self.verticalLayout = QVBoxLayout(self)
        self.controlpanelLayout.addLayout(self.environmentLayout,0,0,Qt.AlignTop)
        self.controlpanelLayout.addLayout(self.verticalLayout,0,1,Qt.AlignCenter)
        self.controlpanelLayout.setColumnStretch(0, 1)
        self.controlpanelLayout.setColumnStretch(1, 3)
        
        self.l1 =  QVBoxLayout()
        self.l1.setAlignment(Qt.AlignTop)

        self.sub1 = QHBoxLayout()
        self.sub1.setAlignment(Qt.AlignLeft)
        self.CH4_txt = QLabel("CH4 (uM)", self)

        self.CH4_box = QtWidgets.QDoubleSpinBox(self)
        self.CH4_box.setSingleStep(0.01)
        self.CH4_box.setDecimals(2)

        self.sub1.addWidget(self.CH4_txt)
        self.sub1.addWidget(self.CH4_box)

        self.l1.addLayout(self.sub1)

        self.CH4 = QSlider(self)
        self.CH4.setOrientation(Qt.Horizontal)
        self.l1.addWidget(self.CH4)

          
        self.l2 =  QVBoxLayout()
        self.l2.setAlignment(Qt.AlignTop)

        self.sub2 = QHBoxLayout()
        self.sub2.setAlignment(Qt.AlignLeft)
        self.H2_txt = QLabel("H2 (mM)", self)

        self.H2_box = QtWidgets.QDoubleSpinBox(self)
        self.H2_box.setSingleStep(0.01)
        self.H2_box.setDecimals(2)

        self.sub2.addWidget(self.H2_txt)
        self.sub2.addWidget(self.H2_box)

        self.l2.addLayout(self.sub2)

        self.H2 = QSlider(self)
        self.H2.setOrientation(Qt.Horizontal)
        self.l2.addWidget(self.H2)


        
        self.l3 =  QVBoxLayout()
        self.l3.setAlignment(Qt.AlignTop)

        self.sub3 = QHBoxLayout()
        self.sub3.setAlignment(Qt.AlignLeft)
        self.DIC_txt = QLabel("DIC (mM)", self)

        self.DIC_box = QtWidgets.QDoubleSpinBox(self)
        self.DIC_box.setSingleStep(0.01)
        self.DIC_box.setDecimals(2)

        self.sub3.addWidget(self.DIC_txt)
        self.sub3.addWidget(self.DIC_box)

        self.l3.addLayout(self.sub3)

        self.DIC = QSlider(self)
        self.DIC.setOrientation(Qt.Horizontal)
        self.l3.addWidget(self.DIC)


        
                
        self.l4 =  QVBoxLayout()
        self.l4.setAlignment(Qt.AlignTop)

        self.sub4 = QHBoxLayout()
        self.sub4.setAlignment(Qt.AlignLeft)
        self.pH_txt = QLabel("pH", self)

        self.pH_box = QtWidgets.QDoubleSpinBox(self)
        self.pH_box.setSingleStep(0.01)
        self.pH_box.setDecimals(2)

        self.sub4.addWidget(self.pH_txt)
        self.sub4.addWidget(self.pH_box)

        self.l4.addLayout(self.sub4)

        self.pH = QSlider(self)
        self.pH.setOrientation(Qt.Horizontal)
        self.l4.addWidget(self.pH)


                        
        self.l5 =  QVBoxLayout()
        self.l5.setAlignment(Qt.AlignTop)

        self.sub5 = QHBoxLayout()
        self.sub5.setAlignment(Qt.AlignLeft)
        self.temperature_txt = QLabel("Temperature (C)", self)

        self.temperature_box = QtWidgets.QDoubleSpinBox(self)
        self.temperature_box.setSingleStep(0.01)
        self.temperature_box.setDecimals(2)

        self.sub5.addWidget(self.temperature_txt)
        self.sub5.addWidget(self.temperature_box)

        self.l5.addLayout(self.sub5)

        self.temperature = QSlider(self)
        self.temperature.setOrientation(Qt.Horizontal)
        self.l5.addWidget(self.temperature)


        self.l6 =  QVBoxLayout()
        self.l6.setAlignment(Qt.AlignTop)

        self.sub6 = QHBoxLayout()
        self.sub6.setAlignment(Qt.AlignLeft)
        self.y_txt = QLabel("Y-axis", self)

        self.y_box = QtWidgets.QComboBox(self)

        self.y_box.addItems(names)        
        self.sub6.addWidget(self.y_txt)
        self.sub6.addWidget(self.y_box)

        self.l6.addLayout(self.sub6)



        self.l7 =  QVBoxLayout()
        self.l7.setAlignment(Qt.AlignTop)
        self.onlyInt = QtGui.QIntValidator()
        self.sub7 = QHBoxLayout()
        self.sub7.setAlignment(Qt.AlignLeft)
        self.volumesoln_txt = QLabel("Solution Volume (mL)", self)
        
        self.volumesoln_box = QtWidgets.QLineEdit(self)
        self.volumesoln_box.setValidator(self.onlyInt)
        
        self.sub8 = QHBoxLayout()
        self.sub8.setAlignment(Qt.AlignLeft)
        self.volumehead_txt = QLabel("Solution Headspace (mL)", self)
        
        self.volumehead_box = QtWidgets.QLineEdit(self)
        self.volumehead_box.setValidator(self.onlyInt)
        
        self.sub7.addWidget(self.volumesoln_txt)
        self.sub7.addWidget(self.volumesoln_box)
        
        self.sub8.addWidget(self.volumehead_txt)
        self.sub8.addWidget(self.volumehead_box)

        self.l7.addLayout(self.sub7)
        self.l7.addLayout(self.sub8)


        #Right Panel
        self.verticalLayout.addLayout(self.l1)
        self.verticalLayout.addLayout(self.l2)
        self.verticalLayout.addLayout(self.l3)
        self.verticalLayout.addLayout(self.l4)

        
        #Left Panel
        self.environmentLayout.addLayout(self.l7)
        self.environmentLayout.addLayout(self.l5)
        self.environmentLayout.addLayout(self.l6)




        self.CH4.valueChanged.connect(lambda: self.setValues('CH4_slider'))
        self.CH4_box.valueChanged.connect(lambda: self.setValues('CH4_box'))

        self.H2.valueChanged.connect(lambda: self.setValues('H2_slider'))
        self.H2_box.valueChanged.connect(lambda: self.setValues('H2_box'))
        
        self.DIC.valueChanged.connect(lambda: self.setValues('DIC_slider'))
        self.DIC_box.valueChanged.connect(lambda: self.setValues('DIC_box'))
        
        self.pH.valueChanged.connect(lambda: self.setValues('pH_slider'))
        self.pH_box.valueChanged.connect(lambda: self.setValues('pH_box'))
        
        self.temperature.valueChanged.connect(lambda: self.setValues('temperature_slider'))
        self.temperature_box.valueChanged.connect(lambda: self.setValues('temperature_box'))
        
        self.volumesoln_box.textChanged.connect(lambda: self.setValues('volumesoln_box'))
        
        self.volumehead_box.textChanged.connect(lambda: self.setValues('volumehead_box'))
        
        # STYLING
        self.CH4_txt.setStyleSheet(QLabel_style)
        self.H2_txt.setStyleSheet(QLabel_style)
        self.DIC_txt.setStyleSheet(QLabel_style)
        self.pH_txt.setStyleSheet(QLabel_style)
        self.temperature_txt.setStyleSheet(QLabel_style)
        self.y_txt.setStyleSheet(QLabel_style)
        self.volumehead_txt.setStyleSheet(QLabel_style)
        self.volumesoln_txt.setStyleSheet(QLabel_style)

    def setValues(self, kind=''):
        if kind == 'CH4_slider':
            self.CH4val = linear_interp(self.CH4.value(),
                                            0, 99, 0.01, 10, lim=1, dec=2)
            self.CH4_box.blockSignals(True)
            self.CH4_box.setValue(self.CH4val)
            self.CH4_box.blockSignals(False)

        elif kind == 'CH4_box':
            self.CH4val = linear_interp(self.CH4_box.value(),
                                            0.01, 10, 0.01, 10, lim=1, dec=2)
            self.CH4.blockSignals(True)
            self.CH4.setValue( linear_interp(self.CH4val,
                                            0.01, 10, 0, 99, lim=1, dec=2) )
            self.CH4.blockSignals(False)
            
        elif kind == 'H2_slider':
            self.H2val = linear_interp(self.H2.value(),
                                            0, 99, 0.01, 10, lim=1, dec=2)
            self.H2_box.blockSignals(True)
            self.H2_box.setValue(self.H2val)
            self.H2_box.blockSignals(False)

        elif kind == 'H2_box':
            self.H2val = linear_interp(self.H2_box.value(),
                                            0.01, 10, 0.01, 10, lim=1, dec=2)
            self.H2.blockSignals(True)
            self.H2.setValue( linear_interp(self.H2val,
                                            0.01, 10, 0, 99, lim=1, dec=2) )
            self.H2.blockSignals(False)
            
        elif kind == 'DIC_slider':
            self.DICval = linear_interp(self.DIC.value(),
                                            0, 99, 0.01, 10, lim=1, dec=2)
            self.DIC_box.blockSignals(True)
            self.DIC_box.setValue(self.DICval)
            self.DIC_box.blockSignals(False)

        elif kind == 'DIC_box':
            self.DICval = linear_interp(self.DIC_box.value(),
                                            0.01, 10, 0.01, 10, lim=1, dec=2)
            self.DIC.blockSignals(True)
            self.DIC.setValue( linear_interp(self.DICval,
                                            0.01, 10, 0, 99, lim=1, dec=2) )
            self.DIC.blockSignals(False)
            
        elif kind == 'pH_slider':
            self.pHval = linear_interp(self.pH.value(),
                                            0, 99, 0.01, 13.99, lim=1, dec=2)
            self.pH_box.blockSignals(True)
            self.pH_box.setValue(self.pHval)
            self.pH_box.blockSignals(False)

        elif kind == 'pH_box':
            self.pHval = linear_interp(self.pH_box.value(),
                                            0.01, 13.99, 0, 13.99, lim=1, dec=2)
            self.pH.blockSignals(True)
            self.pH.setValue( linear_interp(self.pHval,
                                            0.01, 13.99, 0, 99, lim=1, dec=2) )
            self.pH.blockSignals(False)
        elif kind == 'temperature_slider':
            self.temperatureval = linear_interp(self.temperature.value(),
                                            0.01, 99, 0, 100, lim=1, dec=2)
            self.temperature_box.blockSignals(True)
            self.temperature_box.setValue(self.temperatureval)
            self.temperature_box.blockSignals(False)

        elif kind == 'temperature_box':
            self.temperatureval = linear_interp(self.temperature_box.value(),
                                            0, 100, 0, 100, lim=1, dec=2)
            self.temperature.blockSignals(True)
            self.temperature.setValue( linear_interp(self.temperatureval,
                                            0, 100, 0, 99, lim=1, dec=2) )
            self.temperature.blockSignals(False)
            
        elif kind == 'volumesoln_box':
            self.volumesolnval = float(self.volumesoln_box.text())
        
        elif kind == 'volumehead_box':
            self.volumeheadval = float(self.volumehead_box.text())

        else:
            # defaults
            self.CH4_box.setValue(1.0)
            self.H2_box.setValue(1.0)
            self.DIC_box.setValue(2.0)
            self.pH_box.setValue(7.5)
            self.temperature_box.setValue(40)
            self.volumesoln_box.setText("80")
            self.volumehead_box.setText("20")
            self.setValues('CH4_box')
            self.setValues('H2_box')
            self.setValues('DIC_box')
            self.setValues('pH_box')
            self.setValues('temperature_box')
            self.setValues('volumesoln_box')
            self.setValues('volumehead_box')
            
class Widget(QWidget):
    def __init__(self, app, parent=None):
        super(Widget, self).__init__(parent=parent)
                
        self.setStyleSheet(f"Widget {{ background-color: {colors['dark']}; }}")
        
        self.app = app

        self.horizontalLayout = QHBoxLayout(self)
        self.controls = Controls(parent=self)

        self.controls.setValues()
        self.horizontalLayout.addWidget(self.controls)

        self.win = pg.GraphicsLayoutWidget()
        self.horizontalLayout.addWidget(self.win)
        self.plots = [
                        self.win.addPlot(col=1, title="placeholder",
                                        labels={'left':"y",
                                                'bottom':"x"})
                        ]
        
                
        self.update_plot()
        #Set up a timer to refresh plot every 1 ms
        #Implemented to prevent crashing when plots updated too quickly
        self.timer = QtCore.QTimer()
        self.timer.setInterval(1)
        self.timer.timeout.connect(self.update_plot)
        self.timer.start()

    def clear(self):
        print('Clear')
        for p in self.plots:
            p.clear()
        self.update_plot()

    def update_plot(self):

            self.redraw_plots()
            self.win.show()


    def redraw_plots(self):

        CH4_val = self.controls.CH4val
        H2_val = self.controls.H2val
        DIC_val = self.controls.DICval
        pH_val = self.controls.pHval
        temperature_val = self.controls.temperatureval
        y_val = self.controls.y_box.currentText()
        volumesoln_val = self.controls.volumesolnval
        volumehead_val = self.controls.volumeheadval
        
        x_vals, y_vals = get_vals(fac,y=y_val,CH4=CH4_val,H2=H2_val,DIC=DIC_val,
                                  pH=pH_val,temperature=temperature_val,
                                  volumesoln=volumesoln_val,volumehead=volumehead_val)
        
        scatter(self.plots[0],x_vals,y_vals)
        self.plots[0].setTitle(y_val)
        self.plots[0].setVisible(True)
    
    #Timer still running/saved in eiter app or w and causes the script to error
    #out if the script is run again after stopping. Currently a work-around, might
    #want to find a better solution
    def closeEvent(self, evt):
        print("Ended")
        self.timer.stop()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = Widget(app)
    w.show()
    sys.exit(app.exec_())