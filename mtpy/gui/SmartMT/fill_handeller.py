"""
    Description:
        todo to be written

    Usage:
        todo tobe written

    Author: YingzhiGou
    Date: 20/06/2017
"""

import os

import mtpy.core.mt as mt

DEFAULT_GROUP_PROFIX = 'Group'


class FileHandeller:
    """
        Description:
            container that holds all file references and MT object created from the files
    """

    def __init__(self):
        self._file_dict = dict()
        self._group_dict = dict()
        pass

    def add_file(self, file_name, group_id=None):
        """

        :param file_name:
        :param group_id:
        :return:
        """
        file_ref, mt_obj = None
        if isinstance(file_name, str):
            if os.path.isfile(file_name):
                if file_name in self._file_dict and self._file_dict[file_name] is not None:
                    raise FileHandlingException("File %s already loaded." % file_name)
                else:
                    file_ref = file_name
                    mt_obj = mt.MT(file)
        elif isinstance(file_name, mt.MT):
            mt_obj = file_name
            file_ref = file_name.fn
        else:
            raise FileHandlingException("Unsupported input type %s" % type(file_name))

        # add file in to container
        self._file_dict[file_ref] = mt_obj
        # add file to group
        if group_id:
            if isinstance(group_id, list):
                for gid in group_id:
                    self.add_to_group(gid, file_ref)
            else:
                self.add_to_group(group_id, file_ref)
        return True

    def add_files(self, file_list, group_id=None):
        for file_name in file_list:
            self.add_file(file_name, group_id)
        return True

    def add_to_group(self, group_id, file_ref):
        """
        add a file ref to a group
        :param group_id:
        :type group_id str
        :param file_ref:
        :type file_ref str
        :return: True
        """
        if group_id not in self._group_dict:
            self._group_dict[group_id] = set()
        self._group_dict[group_id].update(file_ref)
        return True

    def remove_file_from_group(self, group_id, file_ref):
        """

        :param group_id:
        :type group_id str
        :param file_ref:
        :type file_ref str
        :return:
        """
        if group_id in self._group_dict:
            try:
                self._group_dict[group_id].remove(file_ref)
                return True
            except KeyError:
                return False
        return False


class FileHandlingException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(*args, **kwargs)
