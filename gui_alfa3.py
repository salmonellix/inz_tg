# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\hp\Desktop\inż\myPlug2\gui_alfa3.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(671, 661)
        Form.setLayoutDirection(QtCore.Qt.LeftToRight)
        Form.setStyleSheet("background-color: rgb(44, 44, 44);")
        self.title = QtWidgets.QLabel(Form)
        self.title.setGeometry(QtCore.QRect(230, 10, 251, 51))
        self.title.setStyleSheet("font: 75 16pt \"Myanmar Text\";\n"
"color: rgb(255, 255, 255);\n"
"")
        self.title.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.title.setLineWidth(1)
        self.title.setMidLineWidth(0)
        self.title.setObjectName("title")
        self.line = QtWidgets.QFrame(Form)
        self.line.setGeometry(QtCore.QRect(20, 320, 291, 21))
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.line_2 = QtWidgets.QFrame(Form)
        self.line_2.setGeometry(QtCore.QRect(20, 410, 281, 20))
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.step3 = QtWidgets.QLabel(Form)
        self.step3.setGeometry(QtCore.QRect(90, 440, 211, 16))
        self.step3.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 11pt \"MS Shell Dlg 2\";")
        self.step3.setObjectName("step3")
        self.browse_button1 = QtWidgets.QPushButton(Form)
        self.browse_button1.setGeometry(QtCore.QRect(10, 120, 51, 51))
        self.browse_button1.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255)")
        self.browse_button1.setIconSize(QtCore.QSize(33, 33))
        self.browse_button1.setObjectName("browse_button1")
        self.step2 = QtWidgets.QLabel(Form)
        self.step2.setGeometry(QtCore.QRect(90, 340, 221, 18))
        self.step2.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 11pt \"MS Shell Dlg 2\";")
        self.step2.setObjectName("step2")
        self.step6 = QtWidgets.QLabel(Form)
        self.step6.setGeometry(QtCore.QRect(400, 470, 221, 18))
        self.step6.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 11pt \"MS Shell Dlg 2\";")
        self.step6.setObjectName("step6")
        self.line_3 = QtWidgets.QFrame(Form)
        self.line_3.setGeometry(QtCore.QRect(370, 440, 281, 20))
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.horizontalLayoutWidget_3 = QtWidgets.QWidget(Form)
        self.horizontalLayoutWidget_3.setGeometry(QtCore.QRect(370, 200, 261, 31))
        self.horizontalLayoutWidget_3.setObjectName("horizontalLayoutWidget_3")
        self.box_step4_1 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_3)
        self.box_step4_1.setContentsMargins(0, 0, 0, 0)
        self.box_step4_1.setObjectName("box_step4_1")
        self.label_radius1 = QtWidgets.QLabel(self.horizontalLayoutWidget_3)
        self.label_radius1.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 7pt \"MS Shell Dlg 2\";\n"
"")
        self.label_radius1.setObjectName("label_radius1")
        self.box_step4_1.addWidget(self.label_radius1)
        self.radius1 = QtWidgets.QDoubleSpinBox(self.horizontalLayoutWidget_3)
        self.radius1.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.radius1.sizePolicy().hasHeightForWidth())
        self.radius1.setSizePolicy(sizePolicy)
        self.radius1.setMouseTracking(True)
        self.radius1.setAcceptDrops(True)
        self.radius1.setStyleSheet("background-color: rgb(208, 208, 208);\n"
"selection-color: rgb(90, 181, 249);")
        self.radius1.setAccelerated(True)
        self.radius1.setProperty("value", 3.0)
        self.radius1.setObjectName("radius1")
        self.box_step4_1.addWidget(self.radius1)
        self.ok1 = QtWidgets.QPushButton(self.horizontalLayoutWidget_3)
        self.ok1.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255);")
        self.ok1.setObjectName("ok1")
        self.box_step4_1.addWidget(self.ok1)
        self.horizontalLayoutWidget_4 = QtWidgets.QWidget(Form)
        self.horizontalLayoutWidget_4.setGeometry(QtCore.QRect(370, 500, 261, 31))
        self.horizontalLayoutWidget_4.setObjectName("horizontalLayoutWidget_4")
        self.box_step6 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_4)
        self.box_step6.setContentsMargins(0, 0, 0, 0)
        self.box_step6.setObjectName("box_step6")
        self.label_radius2 = QtWidgets.QLabel(self.horizontalLayoutWidget_4)
        self.label_radius2.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 7pt \"MS Shell Dlg 2\";\n"
"")
        self.label_radius2.setObjectName("label_radius2")
        self.box_step6.addWidget(self.label_radius2)
        self.radius2 = QtWidgets.QDoubleSpinBox(self.horizontalLayoutWidget_4)
        self.radius2.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.radius2.sizePolicy().hasHeightForWidth())
        self.radius2.setSizePolicy(sizePolicy)
        self.radius2.setMouseTracking(True)
        self.radius2.setAcceptDrops(True)
        self.radius2.setStyleSheet("background-color: rgb(208, 208, 208);\n"
"selection-color: rgb(90, 181, 249);")
        self.radius2.setAccelerated(True)
        self.radius2.setProperty("value", 3.0)
        self.radius2.setObjectName("radius2")
        self.box_step6.addWidget(self.radius2)
        self.ok2 = QtWidgets.QPushButton(self.horizontalLayoutWidget_4)
        self.ok2.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255);")
        self.ok2.setObjectName("ok2")
        self.box_step6.addWidget(self.ok2)
        self.button_close = QtWidgets.QPushButton(Form)
        self.button_close.setGeometry(QtCore.QRect(430, 620, 121, 31))
        self.button_close.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255);")
        self.button_close.setObjectName("button_close")
        self.add_hs = QtWidgets.QPushButton(Form)
        self.add_hs.setGeometry(QtCore.QRect(480, 390, 71, 31))
        self.add_hs.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255)")
        self.add_hs.setIconSize(QtCore.QSize(33, 30))
        self.add_hs.setObjectName("add_hs")
        self.line_5 = QtWidgets.QFrame(Form)
        self.line_5.setGeometry(QtCore.QRect(30, 590, 611, 20))
        self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_5.setObjectName("line_5")
        self.step1 = QtWidgets.QLabel(Form)
        self.step1.setGeometry(QtCore.QRect(100, 70, 219, 21))
        self.step1.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 11pt \"MS Shell Dlg 2\";")
        self.step1.setObjectName("step1")
        self.step4 = QtWidgets.QLabel(Form)
        self.step4.setGeometry(QtCore.QRect(400, 70, 219, 21))
        self.step4.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 11pt \"MS Shell Dlg 2\";")
        self.step4.setObjectName("step4")
        self.line_6 = QtWidgets.QFrame(Form)
        self.line_6.setGeometry(QtCore.QRect(323, 100, 20, 461))
        self.line_6.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_6.setObjectName("line_6")
        self.line_7 = QtWidgets.QFrame(Form)
        self.line_7.setGeometry(QtCore.QRect(20, 50, 611, 20))
        self.line_7.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_7.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_7.setObjectName("line_7")
        self.line_8 = QtWidgets.QFrame(Form)
        self.line_8.setGeometry(QtCore.QRect(370, 320, 281, 21))
        self.line_8.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_8.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_8.setObjectName("line_8")
        self.step5 = QtWidgets.QLabel(Form)
        self.step5.setGeometry(QtCore.QRect(400, 360, 221, 18))
        self.step5.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 11pt \"MS Shell Dlg 2\";")
        self.step5.setObjectName("step5")
        self.browser_files1 = QtWidgets.QTextBrowser(Form)
        self.browser_files1.setGeometry(QtCore.QRect(70, 120, 101, 111))
        self.browser_files1.setStyleSheet("background-color: rgb(208, 208, 208);\n"
"selection-color: rgb(118, 197, 230);")
        self.browser_files1.setObjectName("browser_files1")
        self.load1 = QtWidgets.QPushButton(Form)
        self.load1.setGeometry(QtCore.QRect(70, 300, 241, 23))
        self.load1.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255);")
        self.load1.setObjectName("load1")
        self.browser_files1_2 = QtWidgets.QTextBrowser(Form)
        self.browser_files1_2.setGeometry(QtCore.QRect(210, 120, 101, 111))
        self.browser_files1_2.setStyleSheet("background-color: rgb(208, 208, 208);\n"
"selection-color: rgb(118, 197, 230);")
        self.browser_files1_2.setObjectName("browser_files1_2")
        self.label = QtWidgets.QLabel(Form)
        self.label.setGeometry(QtCore.QRect(70, 100, 101, 20))
        self.label.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 8pt \"MS Shell Dlg 2\";")
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(Form)
        self.label_2.setGeometry(QtCore.QRect(210, 100, 101, 20))
        self.label_2.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 8pt \"MS Shell Dlg 2\";")
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(Form)
        self.label_3.setGeometry(QtCore.QRect(180, 160, 21, 21))
        self.label_3.setStyleSheet("font: 75 10pt \"MS Shell Dlg 2\";\n"
"color: rgb(255, 255, 255);")
        self.label_3.setObjectName("label_3")
        self.clear_files1 = QtWidgets.QPushButton(Form)
        self.clear_files1.setGeometry(QtCore.QRect(10, 180, 51, 51))
        self.clear_files1.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255)")
        self.clear_files1.setIconSize(QtCore.QSize(33, 33))
        self.clear_files1.setObjectName("clear_files1")
        self.color_button = QtWidgets.QPushButton(Form)
        self.color_button.setGeometry(QtCore.QRect(120, 550, 121, 31))
        self.color_button.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255);")
        self.color_button.setObjectName("color_button")
        self.clear_all = QtWidgets.QPushButton(Form)
        self.clear_all.setGeometry(QtCore.QRect(120, 620, 121, 31))
        self.clear_all.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255);")
        self.clear_all.setObjectName("clear_all")
        self.hs_browse = QtWidgets.QTextBrowser(Form)
        self.hs_browse.setGeometry(QtCore.QRect(370, 390, 101, 31))
        self.hs_browse.setStyleSheet("background-color: rgb(208, 208, 208);\n"
"selection-color: rgb(118, 197, 230);")
        self.hs_browse.setObjectName("hs_browse")
        self.clear_hs = QtWidgets.QPushButton(Form)
        self.clear_hs.setGeometry(QtCore.QRect(560, 390, 71, 31))
        self.clear_hs.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255)")
        self.clear_hs.setIconSize(QtCore.QSize(33, 30))
        self.clear_hs.setObjectName("clear_hs")
        self.save_file1 = QtWidgets.QPushButton(Form)
        self.save_file1.setGeometry(QtCore.QRect(120, 370, 121, 31))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.save_file1.sizePolicy().hasHeightForWidth())
        self.save_file1.setSizePolicy(sizePolicy)
        self.save_file1.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255);")
        self.save_file1.setObjectName("save_file1")
        self.colors = QtWidgets.QListWidget(Form)
        self.colors.setGeometry(QtCore.QRect(70, 470, 241, 71))
        self.colors.setStyleSheet("background-color: rgb(208, 208, 208);\n"
"color: rgb(255, 255, 255);\n"
"background-color: rgb(44, 44, 44);\n"
"selection-color: rgb(118, 197, 230);")
        self.colors.setObjectName("colors")
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.colors.addItem(item)
        self.save_file2_2 = QtWidgets.QPushButton(Form)
        self.save_file2_2.setGeometry(QtCore.QRect(430, 260, 121, 31))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.save_file2_2.sizePolicy().hasHeightForWidth())
        self.save_file2_2.setSizePolicy(sizePolicy)
        self.save_file2_2.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255);")
        self.save_file2_2.setObjectName("save_file2_2")
        self.save_file3 = QtWidgets.QPushButton(Form)
        self.save_file3.setGeometry(QtCore.QRect(430, 540, 121, 31))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.save_file3.sizePolicy().hasHeightForWidth())
        self.save_file3.setSizePolicy(sizePolicy)
        self.save_file3.setStyleSheet("background-color: rgb(13, 101, 195);\n"
"color: rgb(255, 255, 255);")
        self.save_file3.setObjectName("save_file3")
        self.groupSolvent = QtWidgets.QGroupBox(Form)
        self.groupSolvent.setGeometry(QtCore.QRect(70, 240, 241, 51))
        self.groupSolvent.setStyleSheet("border-color: rgb(40, 40, 40);")
        self.groupSolvent.setTitle("")
        self.groupSolvent.setObjectName("groupSolvent")
        self.oneSolvent = QtWidgets.QRadioButton(self.groupSolvent)
        self.oneSolvent.setGeometry(QtCore.QRect(10, 20, 101, 16))
        self.oneSolvent.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 8pt \"MS Shell Dlg 2\";")
        self.oneSolvent.setCheckable(True)
        self.oneSolvent.setChecked(False)
        self.oneSolvent.setAutoRepeat(False)
        self.oneSolvent.setObjectName("oneSolvent")
        self.multipleSolvent = QtWidgets.QRadioButton(self.groupSolvent)
        self.multipleSolvent.setGeometry(QtCore.QRect(110, 20, 121, 17))
        self.multipleSolvent.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 8pt \"MS Shell Dlg 2\";")
        self.multipleSolvent.setObjectName("multipleSolvent")
        self.groupSimlification = QtWidgets.QGroupBox(Form)
        self.groupSimlification.setGeometry(QtCore.QRect(370, 120, 261, 51))
        self.groupSimlification.setTitle("")
        self.groupSimlification.setObjectName("groupSimlification")
        self.geometric = QtWidgets.QRadioButton(self.groupSimlification)
        self.geometric.setGeometry(QtCore.QRect(10, 20, 131, 16))
        self.geometric.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 8pt \"MS Shell Dlg 2\";")
        self.geometric.setObjectName("geometric")
        self.weights = QtWidgets.QRadioButton(self.groupSimlification)
        self.weights.setGeometry(QtCore.QRect(160, 20, 82, 17))
        self.weights.setStyleSheet("color: rgb(255, 255, 255);\n"
"font: 75 8pt \"MS Shell Dlg 2\";")
        self.weights.setObjectName("weights")

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.title.setText(_translate("Form", "PyMOL HotSpots Plugin"))
        self.step3.setText(_translate("Form", "STEP 3:  Color Spectrum"))
        self.browse_button1.setText(_translate("Form", "BROWSE"))
        self.step2.setText(_translate("Form", "STEP 2:  Save All HotSpots"))
        self.step6.setText(_translate("Form", "STEP 6:     Search groups"))
        self.label_radius1.setText(_translate("Form", "RADIUS     "))
        self.ok1.setText(_translate("Form", "OK"))
        self.label_radius2.setText(_translate("Form", "RADIUS     "))
        self.ok2.setText(_translate("Form", "OK"))
        self.button_close.setText(_translate("Form", "CLOSE"))
        self.add_hs.setText(_translate("Form", "ADD"))
        self.step1.setText(_translate("Form", "STEP 1:  Select Files"))
        self.step4.setText(_translate("Form", "STEP 4:    HotSpots Simplification"))
        self.step5.setText(_translate("Form", "STEP 5:     Choose HotSpots"))
        self.load1.setText(_translate("Form", "LOAD"))
        self.label.setText(_translate("Form", "        HOTSPOTS"))
        self.label_2.setText(_translate("Form", "        PROTEINS"))
        self.label_3.setText(_translate("Form", "  +"))
        self.clear_files1.setText(_translate("Form", "CLEAR"))
        self.color_button.setText(_translate("Form", "COLOR"))
        self.clear_all.setText(_translate("Form", "CLEAR"))
        self.clear_hs.setText(_translate("Form", "CLEAR"))
        self.save_file1.setText(_translate("Form", "SAVE FILE"))
        __sortingEnabled = self.colors.isSortingEnabled()
        self.colors.setSortingEnabled(False)
        item = self.colors.item(0)
        item.setText(_translate("Form", "wheat_brown"))
        item = self.colors.item(1)
        item.setText(_translate("Form", "lightorange_deepolive"))
        item = self.colors.item(2)
        item.setText(_translate("Form", "sand_chocolate"))
        item = self.colors.item(3)
        item.setText(_translate("Form", "blue_red"))
        item = self.colors.item(4)
        item.setText(_translate("Form", "white_black"))
        item = self.colors.item(5)
        item.setText(_translate("Form", "rainbow"))
        item = self.colors.item(6)
        item.setText(_translate("Form", "slate_density"))
        item = self.colors.item(7)
        item.setText(_translate("Form", "dirtyviolet_violetpurple"))
        item = self.colors.item(8)
        item.setText(_translate("Form", "yellow_orange_red"))
        item = self.colors.item(9)
        item.setText(_translate("Form", "smudge_forest"))
        item = self.colors.item(10)
        item.setText(_translate("Form", "darksalmon_raspberry"))
        item = self.colors.item(11)
        item.setText(_translate("Form", "palecyan_deepteal"))
        item = self.colors.item(12)
        item.setText(_translate("Form", "lightpink_lightmagenta"))
        item = self.colors.item(13)
        item.setText(_translate("Form", "lightorange_orange"))
        item = self.colors.item(14)
        item.setText(_translate("Form", "lightblue_purpleblue"))
        item = self.colors.item(15)
        item.setText(_translate("Form", "limon_splitpea"))
        self.colors.setSortingEnabled(__sortingEnabled)
        self.save_file2_2.setText(_translate("Form", "SAVE FILE"))
        self.save_file3.setText(_translate("Form", "SAVE FILE"))
        self.oneSolvent.setText(_translate("Form", "ONE SOLVENT"))
        self.multipleSolvent.setText(_translate("Form", "MULTIPLE SOLVENT"))
        self.geometric.setText(_translate("Form", "GEOMETRIC CENTRE"))
        self.weights.setText(_translate("Form", "HS WEIGHTS"))
import inz_rc


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())
