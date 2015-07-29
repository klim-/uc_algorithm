import ConfigParser
import os
import ast

from core.Logging import myLogger


class ConfigHandler(object):
    """ This class is responsible for parsing the configuration file in
        order to suggest matrices such as pseudoinverses or
        orthocomplements that are not unique.
    """

    def __init__(self):
        self.config = ConfigParser.ConfigParser()
        self.config.optionxform = str

        __dir__ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        filepath = os.path.abspath(os.path.join(__dir__, 'config/default'))

        try:
            data = self.config.read(filepath)
        except Exception as exc:
            myLogger.warn_message("input error!")
            myLogger.debug_message(str(exc))
            print("Input error")

        # fallback values
        # read config-descriptions from dictionary
        self.descr = {}

        with open('core/config_description.py', 'r') as dict:
            data = dict.read()
            self.descr = ast.literal_eval(data)


    def cancle_and_reload(self):
        self.__init__()


    #~ def write(self, section, variable, new_value):
        #~ """ this function saves a variable to config
        #~ """
        #~ self.config.set(str(section), str(variable), str(new_value))

    def read(self, section, variable):
        if self.config.has_option(str(section), str(variable)):
            value = self.config.get(str(section), str(variable))
            # TODO: check if the config value is of same type as fallback value
            myLogger.debug_message(str(variable) + "\": " + str(value) + " (config)")
            return value
        elif str(variable) in self.descr:
            # fallback value
            value = self.descr[str(variable)][1]
            myLogger.debug_message(str(variable) + "\": " + str(value) + " (fallback)")
            return value
        else:
            #pass
            myLogger.error_message("Error! A variable was called that does not exist.")

    def get_boolean(self, section, variable):
        if self.config.has_option(str(section), str(variable)):
            value = self.config.getboolean(str(section), str(variable))
            myLogger.debug_message(str(variable) + "\": " + str(value) + " (config)")
            return value
        elif str(variable) in self.decr:
        #self.descr.has_key(str(variable)):
            # fallback value
            value = self.descr[str(variable)][1]
            myLogger.debug_message(str(variable) + "\": " + str(value) + " (fallback)")
            return value
        else:
            #pass
            myLogger.error_message("Error! A variable was called that does not exist.")


    def apply_changes(self):
        # stores temporary config
        with open('config/default', 'wb') as configfile:
            self.config.write(configfile)


# prepare configData for importing
myConfig = ConfigHandler()

