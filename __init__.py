
from __future__ import absolute_import
from __future__ import print_function


# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

import os

from PyQt5.QtWidgets import QFileDialog


def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('hott', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    '''
    Open our custom dialog
    '''
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

    import sys
    # create a new Window
    dialog = QtWidgets.QDialog()
    QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create('Fusion'))
    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'gui_alfa3.ui')
    form = loadUi(uifile, dialog)


    hs_files=[]
    proteines_files=[]
    hs_paths = []
    proteines_paths = []
    jsons_paths=[]


    def browse_filename_hs_protein():
        protein_path = QFileDialog.getOpenFileName(form, "Open File", "C:\\User\\", 'PDB File (*.pdb)')
        hs_path = QFileDialog.getOpenFileName(form, "Open File", "C:\\User\\", 'MOL2 File (*.mol2)')
        hs_path=hs_path[0]
        hs_paths.append(hs_path)
        hs_path=hs_path.replace(os.path.sep, "\\")

        protein_path = protein_path[0]
        proteines_paths.append(protein_path)
        protein_path = protein_path.replace(os.path.sep, "\\")
        jsons_path=QFileDialog.getOpenFileName(form, "Open File", "C:\\User\\", 'json File (*.json)')
        jsons_path = jsons_path[0]
        jsons_paths.append(jsons_path)

        if hs_path:
            #cmd.load(hs_path)
            hs_name = os.path.basename(os.path.normpath(hs_path))
            hs_name = os.path.splitext(hs_name)[0]
            hs_files.append(hs_name)
            form.browser_files1.setText('\n '.join(hs_files))

            #cmd.load(protein_path)
            protein_name = os.path.basename(os.path.normpath(protein_path))
            protein_name = os.path.splitext(protein_name)[0]
            proteines_files.append(protein_name)
            form.browser_files1_2.setText('\n '.join(proteines_files))




    #load files that are browsed
    def load_hs_protein():
        if len(hs_files) == len(proteines_files):
            cmd.delete('all')
            for nbr in range(0,len(hs_files)):
                cmd.load(hs_paths[nbr])
                cmd.select(hs_files[nbr])
                rezise(nbr)
                cmd.load(proteines_paths[nbr])
            cmd.extra_fit(','.join(hs_files))
            cmd.extra_fit(','.join(proteines_files))

        else:
            form.browser_files1.setText(', '.join('ERROR -_-'))
            form.browser_files1_2.setText('Number of HotSpots and Structures must be the same!!!')

    # change size of hs depends on reference_density_correction from .json file
    def rezise(nbrs):
        def _sele_exists(sele):
            return sele in cmd.get_names('all')

        def _random_string(n=10):
            return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(n))

        def hs_resize(meta_file, selection):

            # Changes size of hotspots depending on normalized partial_charge values.

            if not _sele_exists(selection):
                raise RuntimeError("Selection \"{}\" does not exists.".format(selection))

            # Find free sele name
            temp_sele = _random_string()
            while _sele_exists(temp_sele):
                temp_sele = _random_string()

            states = cmd.count_states(selection=selection)

            with gzip.open(meta_file) as f:
                ref = float(json.load(f)["reference_density_correction"])

            for state in range(1, states + 1):
                stored.info = []
                cmd.iterate_state(state, selection, "stored.info.append((ID, partial_charge))")

                for id_, partial_charge in stored.info:
                    size = log(partial_charge / ref * 1. + 1)

                    cmd.select(temp_sele, "{} and id {}".format(selection, id_), state=state)
                    cmd.set("sphere_scale", value=size, selection=temp_sele)
                    cmd.alter(temp_sele, "b={}".format(partial_charge))

            cmd.delete(temp_sele)

        hs_resize(jsons_paths[nbrs],hs_files[nbrs])


    # save all selected hot spots
    def save_all_hs():
        cmd.select('sele','resn FIL')
        file_name=QFileDialog.getSaveFileName(form, 'Save File','', 'MOL2 File (*.mol2)')
        file_name=file_name[0]
        cmd.extract('all_FIL','sele')
        cmd.save(file_name, 'all_FIL')



    # group hot spots and show only representation
    def hs_sipmly():
        def _sele_exists(sele):
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

        def find_center(*args, **kwargs):
            radius = form.radius1.value()
            color = kwargs.get("color", "red")

            selections = args

            for selection in selections:
                if not _sele_exists(selection):
                    raise RuntimeError("Selection \"{}\" does not exists.".format(selection))

            states = cmd.count_states(selection=selections[0])

            for state in range(1, states + 1):
                stored.info = []

                for selection in selections:
                    cmd.iterate_state(state, selection,
                                      "stored.info.append([ID, partial_charge, (x,y,z), \"{}\"])".format(selection))

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

                    middle_point = group_center(np.array(coords))
                    cmd.pseudoatom(object="simply",
                                   color=color,
                                   pos=middle_point,
                                   state=state)

                    hs_id = free_hotspot_id(hs, exc_indexes=exc_indexes)

                cmd.show(representation="spheres", selection="simply")

        find_center('all_FIL', color = 'magenta')

    # save simplyfied groups of hot spots
    def save_simply():
        cmd.select('simply')
        file_name=QFileDialog.getSaveFileName(form, 'Save File','', 'MOL2 File (*.mol2)')
        file_name=file_name[0]
        cmd.save(file_name, 'simply')

    # callback for the "Browse" button PDB files
    str_files=[]
    def browse_filename_pdb():
        filepath =QFileDialog.getOpenFileName(form, "Open File", "C:\\User\\",'PDB File (*.pdb)')
        filepath=filepath[0]
        filepath=filepath.replace(os.path.sep, "\\")
        if filepath:
            file_name = os.path.basename(os.path.normpath(filepath))
            file_name = os.path.splitext(file_name)[0]
            str_files.append(file_name)
            cmd.load(filepath)
            print(filepath)

    #align structures
    def align_str():
        print(str_files)
        cmd.extra_fit('\n'.join(str_files))

    #make hs color
    def color_hs():
        try:
            cmd.spectrum("pc", "blue_white_red",'all_FIL')
            #cmd.spectrum("pc", "blue_white_red", 'hotspot'+'*')
        except:
            cmd.spectrum("pc", "blue_white_red", 'hotspot')
            cmd.spectrum("pc", "blue_white_red", 'hs'+'*')

    picked_hs=0
    def choose_hs():
        picked_hs=cmd.count_atoms('sele')
        print(picked_hs)
        #picked_hs.append(hs_names)
        form.hs_browse.setText(str(picked_hs))
        cmd.extract('sele', 'sele')

    def clear_pick_hs():
        picked_hs=0
        form.hs_browse.setText(', '.join(picked_hs))



    def clean_all():
        cmd.delete('all')
        hs_files.clear()
        proteines_files.clear()
        hs_paths.clear()
        proteines_paths.clear()
        jsons_paths.clear()
        form.browser_files1_2.setText(', '.join(proteines_files))
        form.browser_files1.setText(', '.join(proteines_files))





    # hook up button callbacks
    # form.button_ray.clicked.connect(run)
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
    return dialog
