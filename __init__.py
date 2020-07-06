

from __future__ import absolute_import
from __future__ import print_function

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.
from typing import List, Any

from PyQt5.QtWidgets import QFileDialog


def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('HS-Detector', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    global dialog
    if dialog is None:
        dialog = make_dialog()

    dialog.show()


def make_dialog():
    # entry point to PyMOL's API
    from pymol import cmd

    # pymol.Qt provides the PyQt5 interface, but may support PyQt4
    # and/or PySide as well
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi
    from PyQt5 import QtGui
    import os
    from pymol import cmd, stored
    import gzip
    import json
    import random
    import string
    from numpy import log
    import numpy as np
    import copy
    import collections
    # create a new Window
    dialog = QtWidgets.QDialog()
    QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create('Fusion'))
    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'gui_alfa3.ui')
    form = loadUi(uifile, dialog)

    # global variables
    Hs_files = []
    Proteines_files = []
    Hs_paths = []
    Proteines_paths = []
    Jsons_paths = []
    solvent = 'non'

    # load protein with hot spots and json file with reference values
    def browse_filename_hs_protein():
        # load protein
        protein_path = QFileDialog.getOpenFileName(form, "Open File", "C:\\User\\", 'PDB File (*.pdb)')
        protein_path = protein_path[0]
        Proteines_paths.append(protein_path)
        protein_path = protein_path.replace(os.path.sep, "\\")
        # load hs file
        hs_path = QFileDialog.getOpenFileName(form, "Open File", "C:\\User\\", 'MOL2 File (*.mol2)')
        hs_path = hs_path[0]
        Hs_paths.append(hs_path)
        hs_path = hs_path.replace(os.path.sep, "\\")
        # load json
        jsons_path = QFileDialog.getOpenFileName(form, "Open File", "C:\\User\\", 'json File (*.json)')
        jsons_path = jsons_path[0]
        Jsons_paths.append(jsons_path)

        if hs_path:
            hs_name = os.path.basename(os.path.normpath(hs_path))
            hs_name = os.path.splitext(hs_name)[0]
            if hs_name != ' ':
                Hs_files.append(hs_name)
                form.browser_files1.setText('\n '.join(Hs_files))

            protein_name = os.path.basename(os.path.normpath(protein_path))
            protein_name = os.path.splitext(protein_name)[0]
            if protein_name != ' ':
                Proteines_files.append(protein_name)
                form.browser_files1_2.setText('\n '.join(Proteines_files))

    # load files that are browsed
    def load_hs_protein():
        if len(Hs_files) == len(Proteines_files):
            cmd.delete('all')
            for nbr in range(0, len(Hs_files)):
                cmd.load(Hs_paths[nbr])
                cmd.select(Hs_files[nbr])
                cmd.load(Proteines_paths[nbr])
                rezise(nbr)
            try:
                cmd.extra_fit(','.join(Hs_files))
                cmd.extra_fit(','.join(Proteines_files))
            except:
                form.browser_files1.setText('ERROR -_-')
                form.browser_files1_2.setText('Number of HotSpots and Structures must be the same!!!')

        else:
            form.browser_files1.setText('ERROR      -_-')
            form.browser_files1_2.setText('Number of HotSpots and Structures must be the same!!!')

    # change size of hs depends on reference_density_correction from .json file (one fil for one pond calculations)
    def rezise(nbrs):
        def is_sele(sele):
            return sele in cmd.get_names('all')

        def _random_string(n=10):
            return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(n))

        refs = []
        ref = 1
        for i in range(0, len(Jsons_paths)):
            with gzip.open(Jsons_paths[i]) as f:
                refOne = float(json.load(f)["reference_density_correction"])
                refs.append(refOne)

        if form.oneSolvent.isChecked():
            ref = np.median(refs)
        if form.multipleSolvent.isChecked():
            with gzip.open(Jsons_paths[nbrs]) as f:
                ref = float(json.load(f)["reference_density_correction"])
        else:
            ref = np.median(refs)

        def hs_resize(meta_file, selection):

            if not is_sele(selection):
                raise RuntimeError("Selection \"{}\" does not exists.".format(selection))

            # Find free sele name
            temp_sele = _random_string()
            while is_sele(temp_sele):
                temp_sele = _random_string()

            states = cmd.count_states(selection=selection)

            for state in range(1, states + 1):
                stored.info = []
                cmd.iterate_state(state, selection, "stored.info.append((ID, partial_charge))")

                for id_, partial_charge in stored.info:
                    size = log(partial_charge / ref + 1)
                    print(ref)
                    print(size)

                    cmd.select(temp_sele, "{} and id {}".format(selection, id_), state=state)
                    cmd.set("sphere_scale", value=size, selection=temp_sele)
                    cmd.alter(temp_sele, "b={}".format(partial_charge))

            cmd.delete(temp_sele)

        hs_resize(Jsons_paths[nbrs], Hs_files[nbrs])

    # save all hs into one mol2 file
    def save_all_hs():
        cmd.select('sele', 'resn FIL')
        file_name = QFileDialog.getSaveFileName(form, 'Save File', '', 'MOL2 File (*.mol2)')
        file_name = file_name[0]
        cmd.create('all_hs', 'sele')
        cmd.save(file_name, 'all_hs', 0)
        cmd.delete('all_hs')

    # simplification geometric centre / hs weights
    def hs_sipmly():

        cmd.delete('simply')

        def is_sele(sele):
            return sele in cmd.get_names("all")

        def find_all_coords(n, hs_id, exc_indexes=[], coords=[]):

            for node in n[hs_id]:
                if (node[0], node[3]) not in exc_indexes:
                    exc_indexes.append((node[0], node[3]))
                    coords.append(node[2])
                    find_all_coords(n, (node[0], node[3]), exc_indexes=exc_indexes, coords=coords)
            return exc_indexes, coords

        def free_hotspot_id(hs, exc_indexes):

            for h in hs:
                if (h[0], h[3]) not in exc_indexes:
                    return h[0], h[3]

            return None

        def group_center(coords):
            length = coords.shape[0]
            sum_x = np.sum(coords[:, 0])
            sum_y = np.sum(coords[:, 1])
            sum_z = np.sum(coords[:, 2])

            return sum_x / length, sum_y / length, sum_z / length

        def find_center(partial_charges, coords, nbr):
            dict_value_coords = {}
            coords_weight=[]

            for p, c in zip(partial_charges,coords):
                dict_value_coords[p] = c
            od_dict = collections.OrderedDict(sorted(dict_value_coords.items(), reverse=True))
            partial_charges=list(od_dict.keys())
            nbr=len(partial_charges)-nbr

            for i in range(0,nbr):
                try:
                    max_hs = partial_charges[i]
                    x = od_dict[max_hs][0]
                    y = od_dict[max_hs][1]
                    z = od_dict[max_hs][2]
                    coords_weight.append([x, y, z])
                except:
                    pass
            return coords_weight


        def hs_gsimplifier(*args, **kwargs):

            radius = form.radius1.value()
            color = kwargs.get("color", "red")

            selections = args

            for selection in selections:
                if not is_sele(selection):
                    raise RuntimeError("Selection \"{}\" does not exists.".format(selection))

            states = cmd.count_states(selection=selections[0])
            hs_values=[]
            all_coords=[]
            for state in range(1, states + 1):
                stored.info = []

                for selection in selections:
                    cmd.iterate_state(state, selection,
                                      "stored.info.append([ID, partial_charge, (x,y,z), \"{}\"])".format(selection))
                for line in stored.info:
                    partial_charge = line[1]
                    coord = line[2]
                    hs_values.append(partial_charge)
                    all_coords.append(coord)

                hs = np.array(copy.copy(stored.info))

                n = dict()
                for h in hs:
                    dists = np.linalg.norm(np.array(h[2]) - list(zip(*hs))[2], axis=1)
                    # Convert to dict due to same hotspots IDs
                    n[(h[0], h[3])] = hs[dists <= radius]

                exc_indexes = []
                hs_id = free_hotspot_id(hs, exc_indexes=exc_indexes)

                while hs_id is not None:
                    excluded, coords = find_all_coords(n, hs_id, exc_indexes=[], coords=[])
                    exc_indexes.extend(excluded)
                    if form.geometric.isChecked():
                        middle_point = group_center(np.array(coords))
                        cmd.pseudoatom(object="simply",
                                       color=color,
                                       pos=middle_point,
                                       state=state)
                    else:
                        middle_points = find_center(hs_values, all_coords, nbr=len(exc_indexes))
                        for point in middle_points:
                            cmd.pseudoatom(object="simply",
                                           color=color,
                                           pos=point,
                                           state=state)

                    hs_id = free_hotspot_id(hs, exc_indexes=exc_indexes)

                cmd.show(representation="spheres", selection="simply")

        cmd.select('hs2', 'resn FIL')
        hs_gsimplifier('hs2', color='hotpink')

    def save_simply():
        cmd.select('simply')
        file_name = QFileDialog.getSaveFileName(form, 'Save File', '', 'MOL2 File (*.mol2)')
        file_name = file_name[0]
        cmd.save(file_name, 'simply')

    # make hs color
    def color_hs():
        color = form.colors.currentItem()
        color = str(color.text())

        try:
            cmd.spectrum("pc", ('%s' % color), 'sele')

        except:
            cmd.spectrum("pc", ('%s' % color), 'hotspot')
            cmd.spectrum("pc", '%s' % color, 'hs' + '*')



    def choose_hs(picked_hs = 0):
        picked_hs = cmd.count_atoms('sele')
        form.hs_browse.setText(str(picked_hs))
        cmd.create('hs', 'sele')

    def clear_pick_hs():
        picked_hs = 0
        form.hs_browse.setText(str(picked_hs))

    def clean_all():
        cmd.delete('all')
        Hs_files.clear()
        Proteines_files.clear()
        Hs_paths.clear()
        Proteines_paths.clear()
        Jsons_paths.clear()
        form.browser_files1_2.setText(', '.join(Proteines_files))
        form.browser_files1.setText(', '.join(Proteines_files))
        form.hs_browse.setText(' ')

    def solvent_one():
        global solvent
        solvent = 'one'

    def solvent_multiple():
        global solvent
        solvent = 'multiple'

    def search_groups():
        radius2 = form.radius2.value()
        cmd.select('sele_atms', 'polymer.protein within %f of hs' %(radius2))

    def save_atms():
        cmd.select('sele_atms')
        file_name = QFileDialog.getSaveFileName(form, 'Save File', '', 'PDB File (*.pdb)')
        file_name = file_name[0]
        cmd.save(file_name, 'sele_atms', -1)

    # hook up button callbacks
    form.browse_button1.clicked.connect(browse_filename_hs_protein)
    form.load1.clicked.connect(load_hs_protein)
    form.button_close.clicked.connect(dialog.close)
    form.save_file1.clicked.connect(save_all_hs)
    form.color_button.clicked.connect(color_hs)
    form.clear_all.clicked.connect(clean_all)
    form.clear_files1.clicked.connect(clean_all)
    form.ok1.clicked.connect(hs_sipmly)
    form.save_file2_2.clicked.connect(save_simply)
    form.add_hs.clicked.connect(choose_hs)
    form.clear_hs.clicked.connect(clear_pick_hs)
    form.oneSolvent.clicked.connect(solvent_one)
    form.multipleSolvent.clicked.connect(solvent_multiple)
    form.ok2.clicked.connect(search_groups)
    form.save_file3.clicked.connect(save_atms)
    return dialog
