import numpy as np
import sys
import os
import datetime

from PyQt5 import QtCore,QtGui
from PyQt5.QtWidgets import QApplication,QWidget,QPushButton,QGridLayout,QFileDialog,QCheckBox,QMessageBox
import pyqtgraph as pg

import matplotlib.pyplot as plt
from generate_test import gen_test

def generatePgColormap(cm_name):
    pltMap = plt.get_cmap(cm_name)
    colors = [pltMap(i) for i in range(pltMap.N)]
    positions = np.linspace(0, 1, len(colors))
    pgMap = pg.ColorMap(positions, colors)
    return pgMap

def getLutfromCmap(cm_name):
    pgMap = generatePgColormap(cm_name)
    return pgMap.getLookupTable()

class MainAssemblerWindow(QWidget):
    def getCuttTimePrefix(self):
        now = datetime.datetime.now()
        return "_".join(map(str, [now.year, now.month, now.day,
                                      now.hour, now.second]))

    def __init__(self,parent=None):
        super().__init__(parent)
        self.settings = {}
        self.baseDir = ""
        if os.path.exists(sys.argv[0]+".ini"):
            f=open(sys.argv[0]+".ini")
            s=f.readline().strip()
            self.baseDir = s
            f.close()
        self.diapason = None
        self.actions = []
        self.settings = {}
        timePrefix = self.getCuttTimePrefix()
        self.settings["actionsfile"] =  timePrefix + "actionsLog.txt"
        self.settings["mainLogFile"] = timePrefix + "mainLog.txt"
        with open(self.settings["mainLogFile"],"a") as fout:
            fout.write(str(datetime.datetime.now())+" ScaffoldIT started\n")
        self.loadUI()
        self.settings["doNorm"] = self.cbNormToSize.isChecked()
        self.show()

    def loadUI(self):
        #create widgets and connect to functions
        self.bnLoadHiCPro = QPushButton("Load HiCPro Matrix")
        self.bnLoadHiCPro.clicked.connect(self.loadHiCPro)

        self.cbNormToSize = QCheckBox("Normalize to bin size")
        self.cbAutoLog = QCheckBox("AutoLog")
        self.cbAutoLog.setChecked(False)

        self.bnLoadTest =  QPushButton("Load Test data")
        self.bnLoadTest.clicked.connect(self.loadTest)

        self.cbSplitMode = QCheckBox("SplitMode")
        self.cbSplitMode.hide()

        self.cbCross = QCheckBox("Add Cross")
        self.cbCross.setChecked(True)

        self.bnLoadGenome = QPushButton("Load Assmebly")
        self.bnLoadGenome.clicked.connect(self.bnLoadGenomeClick)
        self.bnLoadGenome.hide()

        self.bnSaveGenome = QPushButton("Save Assmebly")
        self.bnSaveGenome.clicked.connect(self.bnSaveGenomeClick)
        self.bnSaveGenome.hide()

        self.bnTest = QPushButton("Test")
        #self.bnTest.clicked.connect(self.test)
        self.bnTest.hide()

        self.plotWidget = pg.PlotWidget(enableMenu=False)
        self.gradientWidget = pg.GradientWidget(orientation="right")
        self.gradientWidget.sigGradientChanged.connect(self.changeGradient)

        #create layout
        self.layout = QGridLayout()

        #customize positions of widgets
        #buttons
        self.layout.addWidget(self.bnLoadHiCPro,1,0)
        self.layout.addWidget(self.bnLoadTest, 1, 1)
        self.layout.addWidget(self.cbNormToSize, 1, 2)
        self.layout.addWidget(self.cbAutoLog, 1, 3)
        self.layout.addWidget(self.bnLoadGenome, 1, 4)
        self.layout.addWidget(self.bnSaveGenome, 1, 5)
        self.layout.addWidget(self.cbSplitMode, 1, 6)
        self.layout.addWidget(self.cbCross, 1, 7)
        self.layout.addWidget(self.bnTest, 1, 8)

        #plot and gradient
        self.layout.addWidget(self.plotWidget,0,0,1,9)
        self.layout.addWidget(self.gradientWidget, 0, 9, 2, 1)
        #add layout
        self.setLayout(self.layout)
        #do some tuning =)
        self.resize(1024,768)
        self.setWindowTitle('ScaffoldIt')

    def update_ini(self,s):
        # Write last dir to file
        if len(s) == 2: #This coud be if s comes from QOpenDialog
            s = s[0]
        f = open(sys.argv[0] + ".ini", "w")
        f.write(os.path.dirname(s))
        f.close()
        self.baseDir = s

    def loadTest(self):
        self.generate_test_data()
        self.init_plot()

    def loadHiCPro(self):
        # Write last dir to file
        self.settings["MatrixFile"] = None
        self.settings["BedFile"] = None
        self.settings["MatrixFile"] = QFileDialog.getOpenFileName(caption="Select Hi-C Pro Matrix file",
                                                                  filter="*.dense",directory=self.baseDir)
        if self.settings["MatrixFile"]:
            #Load matrix
            self.settings["MatrixFile"] = self.settings["MatrixFile"][0]
            self.update_ini(self.settings["MatrixFile"])
            print (self.settings["MatrixFile"])
            self.settings["BedFile"] = QFileDialog.getOpenFileName(caption="Select Hi-C Pro Bed file",
                                                                    directory=os.path.dirname(self.settings["MatrixFile"]),
                                                                    filter="*.bed")
            if self.settings["BedFile"]:
                self.settings["BedFile"] = self.settings["BedFile"][0]
                if self.settings["BedFile"] != None and self.settings["MatrixFile"] != None and \
                    self.settings["BedFile"] != "" and self.settings["MatrixFile"] != "":
                    print (self.settings["BedFile"] , self.settings["MatrixFile"] )
                    self.update_ini(self.settings["BedFile"])
                    self.generate_data_fromHiCPro()
                    self.init_plot()

    def bnSaveGenomeClick(self):
        fname = QFileDialog.getSaveFileName(caption="Save genome",directory=self.baseDir)
        if fname:
            fname = fname[0]
            self.update_ini(fname)
            self.save_genome(fname)

    def bnLoadGenomeClick(self):
        fname = QFileDialog.getOpenFileName(caption="Load genome",directory=self.baseDir)
        if fname:
            fname = fname[0]
            self.load_genome(fname)
            self.update_ini(fname)

    def cbSplitModeClick(self):
        self.diapason = None
        self.deselect_data()

    def generate_test_data(self,data=None,scaffolds=None,scaffoldOrientations=None,scaffoldNames=None):
        if data == None or scaffolds == None:
            print("Generating test data")
            self.data,self.scaffolds,self.scaffoldOrientations,self.scaffoldNames = gen_test()
        else:
            self.data,self.scaffolds,self.scaffoldOrientations,self.scaffoldNames = data,scaffolds,scaffoldOrientations,scaffoldNames

        if sum(self.scaffoldOrientations) != len(self.scaffoldOrientations):
            print ("Changing scaffolds orientation")
            reindex = []
            s=0
            for ind,val in enumerate(self.scaffoldOrientations):
                if val:
                    reindex += list(range(s,s + self.scaffolds[ind]))
                else:
                    reindex += list(range(s+self.scaffolds[ind]-1,s-1,-1))
                s += self.scaffolds[ind]
            temp = np.array(self.data)
            for i in range(len(self.data)):
                for j in range(len(self.data)):
                    self.data[i,j] = temp[reindex[i],reindex[j]]
            del temp


        print("Array allocated")

    def generate_data_fromHiCPro(self):
        self.settings["doNorm"] = self.cbNormToSize.checkState()
        print ("Loading Matrix")
        self.data = np.loadtxt(self.settings["MatrixFile"])
        print("Loading BED file")
        bed = np.genfromtxt(self.settings["BedFile"],
                            dtype=np.dtype([("f0","|U50"),("f1",np.int64),("f2",np.int64)]),
                            #dtype=None,
                            usecols=(0,1,2))
        with open(self.settings["mainLogFile"],"a") as fout:
            fout.write(str(datetime.datetime.now())+"Loading Hi-C Pro data:\n")
            fout.write(self.settings["MatrixFile"]+"\t")
            fout.write(self.settings["BedFile"]+"\n")


        print ("Data loaded, reordering matrix")
        self.scaffolds = []
        self.scaffoldNames = []
        count = 0
        last = bed["f0"][0]
        for scaffold in bed["f0"]:
            if scaffold != last:
                self.scaffoldNames.append(str(last))
                self.scaffolds.append(count)
                count = 1
                last = scaffold
            else:
                count += 1

        self.scaffoldNames.append(str(last))
        self.scaffolds.append(count)

        self.scaffoldOrientations = [1]*len(self.scaffoldNames)
        assert sum(self.scaffolds)==len(self.data)

        if self.settings["doNorm"]:
            coefficients = np.max(bed["f2"]-bed["f1"])/(bed["f2"]-bed["f1"])
            assert np.all(coefficients>0)
            self.data *= coefficients
            self.data *= coefficients[:,None]

        toadd = np.min(self.data[np.nonzero(self.data)])
        self.data = np.log2(self.data+toadd)
        print("Done!")

    def init_plot(self,connect = True):
        #connect - connect plot to click function? - True or Flase
        #if you connect multiple times, function will be executed multiple times
        self.imageItem = pg.ImageItem(self.data)
        self.imageItem.setLookupTable(self.gradientWidget.getLookupTable(254))
        #self.imageItem.setLookupTable(getLutfromCmap("autumn"))
        self.plotWidget.addItem(self.imageItem)
        self.plotWidget.getViewBox().invertY(True)
        self.lines = []
        self.update_borders()
        if connect:
            self.plotWidget.sceneObj.sigMouseClicked.connect(self.onClick)
            self.plotWidget.sceneObj.sigMouseMoved.connect(self.onMove)
            #self.plotWidget.sceneObj.keyPressEvent.connect(self.onKeyPress)

        #now hide "load" buttons
        self.bnLoadTest.setDisabled(True)
        self.bnLoadHiCPro.setDisabled(True)
        self.cbNormToSize.setDisabled(True)
        self.bnLoadGenome.show()
        self.bnSaveGenome.show()
        self.cbSplitMode.show()

    def changeGradient(self):
        try:
            self.imageItem
        except:
            return
        print ("Changing image")
        self.imageItem.setLookupTable(self.gradientWidget.getLookupTable(254))

    def update_borders(self):
        print("Updating borders")
        self.borders = np.cumsum(self.scaffolds)
        print (self.borders)
        print("Removing old items")
        for i in self.lines:
            self.plotWidget.removeItem(i)
        print("Done")
        self.lines=[]
        print("Appending new lines")
        for name in ["left","bottom","right","top"]:
            ax = self.plotWidget.getPlotItem().getAxis(name)
            items =list(zip(self.borders,self.scaffoldNames))
            print (items)
            ax.setTicks([items,[]])
            ax.setZValue(1)
            ax.setGrid(150)
        self.imageItem.setZValue(1)
        #for i in self.borders:
            #self.lines.append(self.plotWidget.addLine(x=i))
            #self.lines.append(self.plotWidget.addLine(y=i))

        #self.plotWidget.getPlotItem().getAxis("bottom").setTicks([
                                                         #[(0,"a"),(10,"b")],
                                                         #[(3,"c"),(4,"d")]
                                                         # ])

        self.plotWidget.getPlotItem().getAxis("bottom").setZValue(1)
        print("Done")

    def plot_cross_line(self,x,y):
        try:
            self.plotWidget.removeItem(self.line_x)
            self.plotWidget.removeItem(self.line_y)
        except:
            pass
        if self.cbCross.isChecked():
            self.line_x = self.plotWidget.addLine(x=x,z=10)
            self.line_y = self.plotWidget.addLine(y=y,z=10)

    def get_selected_scaffold_diapason(self,x):
        t = np.searchsorted(self.borders,x,side="right")
        if t == len(self.borders)-1:
            _diapason = self.borders[-2],self.borders[-1]
        elif t==0:
            _diapason = 0,self.borders[0]
        elif t==len(self.borders):
            _diapason = self.borders[-2], len(self.data)
        else:
            _diapason = self.borders[t-1],self.borders[t]
        return _diapason

    def deselect_data(self):
        try:
            self.plotWidget.removeItem(self.ROI)
        except:
            pass
        self.diapason = None

    def select_data(self,newdiapason):
        print ("Selecting",newdiapason)
        #Remove old selection, if exists
        try:
            self.plotWidget.removeItem(self.ROI)
        except:
            pass

        self.ROI = pg.RectROI(pos=(0,newdiapason[0]),
                              size=(len(self.data),newdiapason[1]-newdiapason[0]),
                              movable=False,
                                pen="w")
        self.plotWidget.addItem(self.ROI)

    def addToSelection(self, newdiapason):
        if self.diapason == None:
            return newdiapason
        else:
            newdiapason2 = min(self.diapason[0],newdiapason[0]),\
                            max(self.diapason[1],newdiapason[1])
            self.select_data(newdiapason2)
            return newdiapason2

    def get_scaffolds_in_diapason(self,diapason):
        if diapason == None or len(diapason)==0:
            print (diapason)
            return -1,-1
        left_boundary = diapason[0]
        if left_boundary == 0:
            from_index = 0
        else:
            from_index = list(self.borders).index(left_boundary)+1
        right_boundary = diapason[1]
        print (self.borders,right_boundary)
        to_index =list(self.borders).index(right_boundary)+1
        assert from_index<=to_index
        return from_index,to_index

    def move_contig(self,d_from,d_to,move_type="upper"):
        self.save_actions()
        def do_reindex(d_from,d_to,maxLen,move_type="exchange"):
            reindex = list(range(0,min(d_from[0],d_to[0])))
            if move_type=="exchange":
                if d_from[0]<d_to[0]:
                    reindex += list(range(d_to[0],d_to[1]))
                    reindex += list(range(d_from[1],d_to[0]))
                    reindex += list(range(d_from[0], d_from[1]))
                else:
                    reindex += list(range(d_from[0],d_from[1]))
                    reindex += list(range(d_to[1],d_from[0]))
                    reindex += list(range(d_to[0], d_to[1]))
                reindex += list(range(max(d_from[1],d_to[1]),maxLen))
            elif move_type=="upper":
                if d_from[0] < d_to[0]:
                    reindex += list(range(d_from[1], d_to[0]))
                    reindex += list(range(d_from[0], d_from[1]))
                    reindex += list(range(d_to[0], maxLen))
                else:
                    reindex += list(range(d_from[0], d_from[1]))
                    reindex += list(range(d_to[0], d_from[0]))
                    reindex += list(range(d_from[1], maxLen))
            else:
                raise
            return reindex

        reindex = do_reindex(d_from,d_to,len(self.data), move_type=move_type)
        assert len(reindex)==len(self.data)

        temp = np.zeros_like(self.data)
        for i in range(len(self.data)):
            for j in range(len(self.data)):
                temp[i,j] = self.data[reindex[i],reindex[j]]

        print ("Moving",d_from,"To",d_to)
        self.data = temp
        self.imageItem.setImage(self.data)

        scaffold_from = self.get_scaffolds_in_diapason(d_from)
        assert scaffold_from[1] - scaffold_from[0] >= 1
        scaffold_to = self.get_scaffolds_in_diapason(d_to)
        assert scaffold_to[1] - scaffold_to[0] == 1
        scaffoldsReindex = do_reindex(scaffold_from,scaffold_to,len(self.scaffolds), move_type=move_type)
        self.actions.append(move_type + "\t" +
                            "+".join([self.scaffoldNames[i] for i in range(scaffold_from[0],scaffold_from[1])]) + \
                            "\t" + \
                            "+".join([self.scaffoldNames[i] for i in range(scaffold_to[0],scaffold_to[1])]))
        assert len(scaffoldsReindex)==len(self.scaffolds)
        self.scaffolds = [self.scaffolds[i] for i in scaffoldsReindex]
        self.scaffoldNames = [self.scaffoldNames[i] for i in scaffoldsReindex]

        self.update_borders()
        self.save_actions()

    def invert_contig(self,diapason):
        self.save_actions()
        print("Inverting contig")
        reindex = list(range(0,diapason[0]))
        reindex += reversed(list(range(diapason[0], diapason[1])))
        reindex += list(range(diapason[1], len(self.data)))
        temp = np.zeros_like(self.data)
        for i in range(len(self.data)):
            for j in range(len(self.data)):
                temp[i,j] = self.data[reindex[i],reindex[j]]

        self.data = temp
        self.imageItem.setImage(self.data)

        scaffold_id = list(self.borders).index(diapason[1])
        if self.scaffoldOrientations[scaffold_id] == 1:
            self.scaffoldOrientations[scaffold_id] = 0
        elif self.scaffoldOrientations[scaffold_id] == 0:
            self.scaffoldOrientations[scaffold_id] = 1
        else:
            raise

        self.actions.append("i\t"+self.scaffoldNames[scaffold_id ])
        self.save_actions()

    def save_actions(self):
        print ("Saving actions...")
        fname = self.settings["actionsfile"]
        with open(fname,"w") as out:
            out.write("\n".join(self.actions)+"\n")
        self.save_genome("genome_log.txt")

    def save_genome(self,fname="genome.txt"):
        def orientation_to_str(i):
            if i==1:
                return  "+"
            elif i==0:
                return  "-"
            else:
                raise
        print ("Saving genome...")

        if self.cbAutoLog.isChecked():
            fname = os.path.join(
                os.path.dirname(fname),
                self.getCuttTimePrefix() + "_" + os.path.basename(fname))

        with open(fname,"w") as out:
            for ind,val in enumerate(self.scaffoldNames):
                out.write(val+"\t"+orientation_to_str(self.scaffoldOrientations[ind])+"\n")
        print ("Data saved to "+fname)
        return fname

    def showErrorMessage(self,msgtext="Unknown Error",msgInfoText="Unknown Error Occur",winTitle="Error"):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText(msgtext)
        msg.setInformativeText(msgInfoText)
        msg.setWindowTitle(winTitle)
        msg.setStandardButtons(QMessageBox.Cancel)
        return msg.exec_()

    def load_genome(self,fname):
        def orientation_to_int(o):
            assert o in ["+","-","1","0"]
            if o == "+":
                return 1
            elif o == "-":
                return 0
            else:
                return int(o)

        try:
            errMsg = "Cannot read assmebly file"
            print (fname)
            new_genome = np.genfromtxt(fname, dtype=np.dtype([("f0", "|U50"), ("f1", "|U50")]))
            errMsg = "Number of scaffold does not match between genomes"
            assert len(new_genome) == len(self.scaffolds)
            errMsg = "Total length of scaffolds does not match matrix"
            assert sum(self.scaffolds) == len(self.data)
            errMsg = "Duplicates found in genome file"
            assert len(np.unique(new_genome["f0"]))==len(new_genome["f0"])

            scaffolds_reindex = []
            new_orientations = []
            reindex = []
            for ind,val in enumerate(new_genome["f0"]):
                if not (val in self.scaffoldNames):
                    errMsg = "Scaffold "+str(val)+" not in scaffolds"
                    raise Exception("GenomeLoadException")
                old_id = self.scaffoldNames.index(val)
                scaffolds_reindex.append(old_id)

                if old_id == 0:
                    scaffold_begin = 0
                else:
                    scaffold_begin = self.borders[old_id-1]
                scaffold_end = self.borders[old_id]

                new_orientations.append(orientation_to_int(new_genome["f1"][ind]))

                if new_orientations[-1] != self.scaffoldOrientations[ind]:
                    reindex += list(range(scaffold_end-1,scaffold_begin-1,-1))
                elif new_orientations[-1] == 1 or new_orientations[-1] == 0:
                    reindex += list(range(scaffold_begin,scaffold_end,1))
                else:
                    print (new_orientations[-1],self.scaffoldOrientations[ind])
                    errMsg = "Unknown orientation for scaffold"+str(val)
                    raise Exception("GenomeLoadException")
            errMsg = "Wrong number/lengths of scaffolds"
            assert len(reindex) == len(self.data)
        except:
            self.showErrorMessage(msgtext="Error loading genome",msgInfoText=errMsg)
            return 0


        temp = np.zeros_like(self.data)
        for i in range(len(self.data)):
            for j in range(len(self.data)):
                temp[i,j] = self.data[reindex[i],reindex[j]]
        self.data = temp
        self.imageItem.setImage(self.data)
        assert len(self.scaffolds) == len(scaffolds_reindex) == len(self.scaffoldNames)
        self.scaffolds = [self.scaffolds[i] for i in scaffolds_reindex]
        self.scaffoldNames = [self.scaffoldNames[i] for i in scaffolds_reindex]
        self.scaffoldOrientations = new_orientations
        self.update_borders()

        with open(self.settings["mainLogFile"],"a") as fout:
            fout.write("Loaded assembly "+fname)

    def split_scaffold(self,scaffold_id,bin):
        try:
            self.settings["BedFile"]
            open(self.settings["BedFile"])
        except:
            self.showErrorMessage("Splitting error","Hi-CPro Bed file not found","Error")
            return

        scaffold_length = self.scaffolds[scaffold_id]
        print ("Scaffold length = ",scaffold_length," bin = ",bin)

        if bin == scaffold_length:
            self.showErrorMessage("Splitting error","Cannot split after last bin","Error")
            return

        if self.scaffoldOrientations[scaffold_id] == 0:
            bin = scaffold_length - bin - 2


        old_scaffolds_N = len(self.scaffoldNames)
        new_scaffold_names = []
        new_scaffold_orient = []
        for ind,(name,orientation) in enumerate(zip(self.scaffoldNames,self.scaffoldOrientations)):
            if ind == scaffold_id:
                if self.scaffoldOrientations[scaffold_id] == 0:
                    new_scaffold_names.append(name + ".2")
                    new_scaffold_names.append(name + ".1")
                else:
                    new_scaffold_names.append(name + ".1")
                    new_scaffold_names.append(name + ".2")
                new_scaffold_orient += [orientation,orientation]
            else:
                new_scaffold_names.append(name)
                new_scaffold_orient.append(orientation)
        assert old_scaffolds_N == len(new_scaffold_names) - 1 == len(new_scaffold_orient) - 1

        with open(self.settings["BedFile"]) as fin,\
            open(self.settings["BedFile"]+".split"+new_scaffold_names[scaffold_id]+".bed","w") as fout:
            curr_ind = 0
            for line in fin:
                if line.startswith("#"):
                    fout.write(line)
                else:
                    data = line.strip().split()
                    if data[0] == self.scaffoldNames[scaffold_id]:
                        if curr_ind <= bin:
                            fout.write(data[0]+".1"+"\t")
                            fout.write("\t".join(data[1:])+"\n")
                        else:
                            fout.write(data[0]+".2"+"\t")
                            fout.write("\t".join(data[1:])+"\n")
                        curr_ind += 1
                    else:
                        fout.write(line)
            assert curr_ind == self.scaffolds[scaffold_id]

        self.scaffoldOrientations = new_scaffold_orient
        self.scaffoldNames = new_scaffold_names
        genomefile = self.save_genome(fname=self.settings["BedFile"]+"split"+new_scaffold_names[scaffold_id]+".genome")

        self.settings["BedFile"] = self.settings["BedFile"]+".split"+new_scaffold_names[scaffold_id]+".bed"
        self.generate_data_fromHiCPro()
        self.plotWidget.clear()
        self.init_plot(connect=False)
        self.load_genome(genomefile)

    def onMove(self,point):
        #posInView = self.plotWidget.getViewBox().mapSceneToView(event.scenePos())
        print (point.x(),point.y())
        posInView = self.plotWidget.getViewBox().mapSceneToView(point)
        print (posInView.x(),posInView.y())
        self.Pos = posInView
        self.plot_cross_line(posInView.x(),posInView.y())

    def onClick(self,event):
        print("--------")
        modifiers = QtGui.QApplication.keyboardModifiers()
        #print (event.button())
        posInView = self.plotWidget.getViewBox().mapSceneToView(event.scenePos())
        newdiapason = self.get_selected_scaffold_diapason(int(posInView.y()))

        if self.cbSplitMode.isChecked(): #click done in splitMode
            sel_scaffolds = self.get_scaffolds_in_diapason(newdiapason)
            if sel_scaffolds[1] - sel_scaffolds[0] != 1:
                self.showErrorMessage("Splitting error","Wrong selection","Splitting error")
                return
            print ("ScaffNames=",self.scaffoldNames)
            print (sel_scaffolds)
            sel_scaffold_id =  sel_scaffolds[0]
            print (sel_scaffold_id)
            #DEBUG
            print ("Selecting scaffold ",self.scaffoldNames[sel_scaffold_id])
            sel_bin = int(posInView.y())-newdiapason[0]
            if sel_bin >= self.scaffolds[sel_scaffold_id]-1:
                self.showErrorMessage("Splitting selection error","Cannot select last bin in scaffold","Splitting Error")
                return
            print ("Selecting bin",sel_bin)
            self.split_scaffold(sel_scaffold_id,sel_bin)
            return

        if event.button() == 2:  # right click for inversion
            d = self.get_scaffolds_in_diapason(self.diapason)
            if d[1]-d[0] != 1:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Information)
                msg.setText("Inversion Error")
                msg.setInformativeText("Only single scaffold can be inverted")
                msg.setWindowTitle("Inversion Error")
                msg.setStandardButtons(QMessageBox.Cancel)
                retval = msg.exec_()
                return
            else:
                self.invert_contig(newdiapason)
                self.diapason = None
                self.deselect_data()
                return

        if modifiers == QtCore.Qt.ShiftModifier:
            print ("Adding")
            self.diapason = self.addToSelection(newdiapason)
            #DEBUG
            print(self.diapason)
        else:
            self.select_data(newdiapason)
            if (self.diapason != None) and (newdiapason != self.diapason) and \
                    not (newdiapason[0]>=self.diapason[0] and newdiapason[1]<=self.diapason[1]):
                print ("Moving contigs",self.diapason,newdiapason)
                self.move_contig(self.diapason,newdiapason)
                self.diapason = None
                self.deselect_data()
            elif self.diapason == None or \
                    (newdiapason[0] > self.diapason[0] and newdiapason[1] < self.diapason[1]):
                print ("Selecting diapason")
                self.diapason = newdiapason
                self.select_data(newdiapason)
            elif self.diapason == newdiapason:
                print ("DeSelecting diapason")
                self.deselect_data()
                self.diapason = None
        print ("Scaffolds in diapason:",self.get_scaffolds_in_diapason(self.diapason))

    def keyPressEvent(self, e):
        if e.key() == QtCore.Qt.Key_Left or e.key() == QtCore.Qt.Key_Right:
            print("Left/Right")
            cursor = self.cursor()
            print(self.Pos.x(), self.Pos.y())
            newPos = self.plotWidget.getViewBox().mapViewToScene(self.Pos)
            print (newPos.x(),newPos.y())
            point  = QtCore.QPoint(newPos.x(),newPos.y())
            print (point)
            print (self.plotWidget.mapToGlobal(point))
            #newPos = self.mapFrom(self.plotWidget,point)
            #print (newPos)
            cursor.setPos(point)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    mainWin = MainAssemblerWindow()
    mainWin.show()
    sys.exit( app.exec_() )