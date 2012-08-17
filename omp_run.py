import datetime
import os
import re
import string
import subprocess
import sys
import time

g_omp_parameters = {}
g_omp_runs = []

def main():
    omp_params_file_path = None
    
    if (len(sys.argv) > 1):
        omp_params_file_path = sys.argv[1]
        
    initRunParameters(omp_params_file_path)    
    runOMP()
    
def initRunParameters(omp_params_file_path=None):
    errors_in_param_file = False
    
    if (omp_params_file_path == None):
        working_dir = sys.path[0].replace("\\", "/")
        omp_params_file_path = working_dir + "/omp_parameters.txt"
    
    if (not os.path.exists(omp_params_file_path)):
        raise StandardError, "The OMP parameters file does not exist."
    
    f = open(omp_params_file_path, 'r')
    datain = f.read()
    f.close()
    
    param_sections_list = datain.split("\n\n")
    for section in param_sections_list:
        section_list = section.split("\n")
        
        if (section_list[0].strip() == "[OMP_Configuation]" or section_list[0].strip() == "[OMP_Experimental_Conditions]"):
            for i in xrange(1,len(section_list)):
                section_item = section_list[i]
                
                if (re.match("^\#", section_item, re.IGNORECASE)):
                    continue
                
                section_item_list = section_item.split("=")
                
                key = section_item_list[0].strip()
                value = section_item_list[1].strip()
                g_omp_parameters[key] = value
        elif (section_list[0].strip() == "[OMP_Run_Info]"):
            omp_run_hash = {}
            for i in xrange(1,len(section_list)):
                section_item = section_list[i]
                
                if (re.match("^\#", section_item, re.IGNORECASE)):
                    continue
                
                section_item_list = section_item.split("=")
                key = section_item_list[0].strip()
                value = section_item_list[1].strip()
                omp_run_hash[key] = value
                
            g_omp_runs.append(omp_run_hash)

    print(len(g_omp_runs))

def runOMP():
    start_time = datetime.datetime.now() 

    exe_file_path = g_omp_parameters["OMP_EXE_FILE_PATH"]
    
    if (g_omp_parameters["OUTPUT_FILE_DELIMITER"] != ""):
        output_file_delimiter = g_omp_parameters["OUTPUT_FILE_DELIMITER"]
    else:
        output_file_delimiter = "\t"
        
    if (g_omp_parameters["INPUT_FILE_DELIMITER"] != ""):
        input_file_delimiter = g_omp_parameters["INPUT_FILE_DELIMITER"]
    else:
        input_file_delimiter = "\t"
    
    for run_info in g_omp_runs:
        seq_input_file_path = run_info["INPUT_FILE_PATH"]
        final_omp_output_file_path = run_info["OUTPUT_FILE_PATH"]
        
        final_omp_output_base_dir = os.path.dirname(final_omp_output_file_path) + "/"
        if (not os.path.exists(final_omp_output_base_dir)):
            os.makedirs(final_omp_output_base_dir)
        
        support_file_dir = run_info["OMP_SUPPORT_FILE_DIR"].strip()
        if (not os.path.exists(support_file_dir)):
            os.makedirs(support_file_dir)
        
        f = open(seq_input_file_path, 'r')
        datain = f.readlines()
        f.close()
        
        dataout = []
        header_list = ["probe_id", "target_id", "species_unit_structure_type", "STRUCTURE_ID", "DELTA-G", 
                       "DELTA-H", "DELTA-S", "MELTING_TEMPERATURE", "CONCENTRATION", "PERCENT_BOUND", 
                        "PERCENT_BOUND(2)", "TARGET_START_POS", "NET_DG", "NET_TM"]
        header = string.join(header_list, output_file_delimiter) + "\n"
        
        f = open(final_omp_output_file_path, 'w')
        f.write(header)
        
        num_rows = len(datain)
        for i in xrange(1,num_rows):
            line = datain[i]
            line = line.strip()
            fields = line.split(input_file_delimiter)
            
            probe_id = fields[0].strip()
            probe_seq = fields[1].strip()
            target_id = fields[2].strip()
            target_seq = fields[3].strip()
            target_type = fields[4].strip()
            save_support_files = fields[5].strip()
            if (save_support_files == "FALSE" or save_support_files == ""):
                save_support_files = False
            elif save_support_files == "TRUE":
                save_support_files = True
            
            if (probe_id == "" or probe_seq == "" or target_id == "" or target_seq == "" or target_type == ""):
                continue
            
            oef_file_name = probe_id + "_vs_" + target_id + ".OEF"
            oef_file_path = final_omp_output_base_dir + oef_file_name
            
            nal_file_name = probe_id + "_vs_" + target_id + ".nal"
            nal_file_path = final_omp_output_base_dir + nal_file_name
            
            oof_file_name = probe_id + "_vs_" + target_id + ".OOF"
            oof_file_path = final_omp_output_base_dir + oof_file_name
            
            tbs_probe_monomer_file_name = probe_id + "_vs_" + target_id + "_" + probe_id + ".tbs"
            tbs_probe_monomer_file_path = final_omp_output_base_dir + tbs_probe_monomer_file_name
            
	    tbs_target_monomer_file_name = probe_id + "_vs_" + target_id + "_" + target_id + ".tbs"
            tbs_target_monomer_file_path = final_omp_output_base_dir + tbs_target_monomer_file_name
            
            tbs_target_homodimer_file_name = probe_id + "_vs_" + target_id + "_" + target_id + "_" + target_id + ".tbs"
            tbs_target_homodimer_file_path = final_omp_output_base_dir + tbs_target_homodimer_file_name
            
            tbs_probe_target_heterodimer_file_name = probe_id + "_vs_" + target_id + "_" + target_id + "_" + probe_id + ".tbs"
            tbs_probe_target_heterodimer_file_path = final_omp_output_base_dir + tbs_probe_target_heterodimer_file_name
	    
	    ta_target_monomer_file_name = probe_id + "_vs_" + target_id + "_" + target_id + ".ta"
            ta_target_monomer_file_path = final_omp_output_base_dir + ta_target_monomer_file_name
            
            createOEF(oef_file_path, probe_id, probe_seq, target_id, target_seq, target_type)
            
            if (g_omp_parameters["TEST_MODE"] == "TRUE"):
                args = [oef_file_path, oof_file_path, probe_id, target_id]
            else:
                args = [oef_file_path, oof_file_path]
                
            try:
                retcode = subprocess.call(exe_file_path + " " + string.join(args, " "), shell=True)
                line = parseOOF(oof_file_path)
                f.write(line)
    #                if retcode < 0:
    #                    print >> sys.stderr, "Child was terminated by signal", -retcode
    #                else:
    #                    print >> sys.stderr, "Child returned", retcode
            except OSError, e:
    #                print >> sys.stderr, "Execution failed:", e
                raise
            finally:
                if (os.path.exists(oef_file_path)):    
                    if (save_support_files and support_file_dir != ""):
                        file_name = os.path.basename(oef_file_path)
                        new_file_path = support_file_dir + file_name
                        os.rename(oef_file_path, new_file_path)
                    else:
                        os.remove(oef_file_path)
                
                if (os.path.exists(nal_file_path)):
                    if (save_support_files and support_file_dir != ""):
                        file_name = os.path.basename(nal_file_path)
                        new_file_path = support_file_dir + file_name
                        os.rename(nal_file_path, new_file_path)
                    else:
                        os.remove(nal_file_path)
                
                if (os.path.exists(oof_file_path)):
                    if (save_support_files and support_file_dir != ""):
                        file_name = os.path.basename(oof_file_path)
                        new_file_path = support_file_dir + file_name
                        os.rename(oof_file_path, new_file_path)
                    else: 
                        os.remove(oof_file_path)
                if (os.path.exists(tbs_probe_monomer_file_path)):
                    if (save_support_files and support_file_dir != ""):
                        file_name = os.path.basename(tbs_probe_monomer_file_path)
                        new_file_name = probe_id + "_" + target_id + "_PROBE_MONOMER.tbs"
                        new_file_path = support_file_dir + new_file_name
                        os.rename(tbs_probe_monomer_file_path, new_file_path)
                    else: 
                        os.remove(tbs_probe_monomer_file_path)
                        
                if (os.path.exists(tbs_target_monomer_file_path)):
		    if (save_support_files and support_file_dir != ""):
                        file_name = os.path.basename(tbs_target_monomer_file_path)
                        new_file_name = probe_id + "_" + target_id + "_TARGET_MONOMER.tbs"
                        new_file_path = support_file_dir + new_file_name
                        os.rename(tbs_target_monomer_file_path, new_file_path)
                    else: 
                        os.remove(tbs_target_monomer_file_path)

		if (os.path.exists(ta_target_monomer_file_path)):
                    if (save_support_files and support_file_dir != ""):
                        file_name = os.path.basename(ta_target_monomer_file_path)
                        new_file_name = probe_id + "_" + target_id + "_TARGET_MONOMER.ta"
                        new_file_path = support_file_dir + new_file_name
                        os.rename(ta_target_monomer_file_path, new_file_path)
                    else: 
                        os.remove(ta_target_monomer_file_path)
                
                if (os.path.exists(tbs_target_homodimer_file_path)):
		    if (save_support_files and support_file_dir != ""):
			file_name = os.path.basename(tbs_target_homodimer_file_path)
                        new_file_name = probe_id + "_" + target_id + "_TARGET_HOMODIMER.tbs"
                        new_file_path = support_file_dir + new_file_name
                        os.rename(tbs_target_homodimer_file_path, new_file_path)
                    else: 
                        os.remove(tbs_target_homodimer_file_path)
                if (os.path.exists(tbs_probe_target_heterodimer_file_path)):
		    if (save_support_files and support_file_dir != ""):
                        file_name = os.path.basename(tbs_probe_target_heterodimer_file_path)
                        new_file_name = probe_id + "_" + target_id + "_PROBE_TARGET_HETERODIMER.tbs"
                        new_file_path = support_file_dir + new_file_name
                        os.rename(tbs_probe_target_heterodimer_file_path, new_file_path)
                    else: 
                        os.remove(tbs_probe_target_heterodimer_file_path)
                    
        end_time = datetime.datetime.now() 
        
        num_probe_target_comparisons = str(num_rows-1)
        time_to_completion = str(end_time-start_time)
        log_line_list = [seq_input_file_path, num_probe_target_comparisons, start_time, end_time, time_to_completion]
        
        f.close()
        
        addLogEntry(log_line_list)
        
        print "Finished OMP run for file '" + seq_input_file_path + "'" 
    
def parseOOF(input_file_path):
    dataout = []
    
    if (g_omp_parameters["OUTPUT_FILE_DELIMITER"] != ""):
        delimiter = g_omp_parameters["OUTPUT_FILE_DELIMITER"]
    else:
        delimiter = "\t"
    
    f = open(input_file_path, 'r')
    oof_data = f.read()
    f.close()
    
    oof_list = oof_data.split("SPECIES=")
    
    structure_fields = ["DELTA-G", "DELTA-H", "DELTA-S", "MELTING_TEMPERATURE", "CONCENTRATION", "PERCENT_BOUND", 
                        "PERCENT_BOUND(2)", "TARGET_START_POS", "NET_DG", "NET_TM"]
    
    for species_unit in oof_list[1:]:
        species_unit = species_unit.strip()
        species_unit_structures = species_unit.split("STRUCTURE=")
        
        species_unit_header = species_unit_structures[0].strip()
        species_unit_header_list = species_unit_header.split("\n")
        species_unit_structure_type = species_unit_header_list[0]
        species_unit_seq_id_list = species_unit_header_list[1].split("=")[1].split("+")
        
        if (species_unit_structure_type == "MONOMER" or species_unit_structure_type == "HOMODIMER"):
            probe_id = species_unit_seq_id_list[0]
            target_id = species_unit_seq_id_list[0]
        elif (species_unit_structure_type == "HETERODIMER"):
            probe_id = species_unit_seq_id_list[1]
            target_id = species_unit_seq_id_list[0]
        
        line_list_base = [probe_id, target_id, species_unit_structure_type]
                
        for structure in species_unit_structures[1:]:
            line_list = []
            line_list.extend(line_list_base)
            structure = structure.split("\n")
            
            structure_value_hash = {}
            structure_value_hash["STRUCTURE_ID"] = structure[0]
            
            for structure_value in structure[1:]:
                structure_value = structure_value.strip()

                if (structure_value != ""):
                    field_name = structure_value.split("=")[0]
                    field_value = structure_value.split("=")[1]
                    structure_value_hash[field_name] = field_value
                    
            line_list.append(structure_value_hash["STRUCTURE_ID"])
            
            for field_name in structure_fields:
                if (structure_value_hash.has_key(field_name)):
                    line_list.append(structure_value_hash[field_name])
                else:
                    line_list.append("\\N")
                    
            line = string.join(line_list, delimiter) + "\n"
            dataout.append(line)
        
    return string.join(dataout, "")
            
def createOEF(file_path, probe_id, probe_seq, target_id, target_seq, target_type):
    dataout = []
    dataout.append("\n[project information]\n")
    dataout.append("GENERATE_LOG=FALSE\n")
    dataout.append("DEFAULT_GENERATE_NETTM=TRUE\n")
    dataout.append("GENERATE_NUMANALY=True\n")
    
    if (g_omp_parameters["NUMANALY_MIN_TEMPERATURE"] != ""):
        dataout.append("NUMANALY_MIN_TEMPERATURE=" + g_omp_parameters["NUMANALY_MIN_TEMPERATURE"] + "\n")
    else:
        dataout.append("NUMANALY_MIN_TEMPERATURE=37\n")
    
    if (g_omp_parameters["NUMANALY_MAX_TEMPERATURE"] != ""):
        dataout.append("NUMANALY_MAX_TEMPERATURE=" + g_omp_parameters["NUMANALY_MAX_TEMPERATURE"] + "\n")
    else:
        dataout.append("NUMANALY_MAX_TEMPERATURE=42\n")
        
    dataout.append("CONVERGENCE_CYCLES=100000\n")
    dataout.append("CONVERGENCE_TIMEOUT=100\n")
    dataout.append("DELTA_CP=0\n")
    dataout.append("CALCULATE_EXTENSIBILITY=TRUE\n")
    dataout.append("EXTENSION_OVERHANG=1\n")
    dataout.append("EXTENSION_WINDOW=3\n")
    dataout.append("MULTIPLEX_ENABLE=FALSE\n")
    dataout.append("MULTIPLEX_ENERGY_THRESHOLD=8\n")
    dataout.append("MULTIPLEX_ENERGY_WINDOW=50\n")
    dataout.append("[solution]\n")
    
    if (g_omp_parameters["ASSAY_TEMPERATURE"] != ""):
        dataout.append("ASSAY_TEMPERATURE=" + g_omp_parameters["ASSAY_TEMPERATURE"] + "\n")
    else:
        dataout.append("ASSAY_TEMPERATURE=37\n")
    
    if (g_omp_parameters["MAGNESIUM_CONCENTRATION"] != ""):
        dataout.append("MAGNESIUM_CONCENTRATION=" + g_omp_parameters["MAGNESIUM_CONCENTRATION"] + "\n")
    else:
        dataout.append("MAGNESIUM_CONCENTRATION=0\n")
    
    if (g_omp_parameters["SODIUM_CONCENTRATION"] != ""):
        dataout.append("SODIUM_CONCENTRATION=" + g_omp_parameters["SODIUM_CONCENTRATION"] + "\n")
    else:
        dataout.append("SODIUM_CONCENTRATION=1\n")
    
    if (g_omp_parameters["GLYCEROL_CONCENTRATION"] != ""):
        dataout.append("GLYCEROL_CONCENTRATION=" + g_omp_parameters["GLYCEROL_CONCENTRATION"] + "\n")
    else:
        dataout.append("GLYCEROL_CONCENTRATION=0\n")
#    dataout.append("GLYCEROL_CONC_UNITS=PERCENT\n") # Not supported
    
    if (g_omp_parameters["DMSO_CONCENTRATION"] != ""):
        dataout.append("DMSO_CONCENTRATION=" + g_omp_parameters["DMSO_CONCENTRATION"] + "\n")
    else:
        dataout.append("DMSO_CONCENTRATION=0\n")
#    dataout.append("DMSO_CONC_UNITS=PERCENT\n") # Not supported
    
    if (g_omp_parameters["FORMAMIDE_CONCENTRATION"] != ""):
        dataout.append("FORMAMIDE_CONCENTRATION=" + g_omp_parameters["FORMAMIDE_CONCENTRATION"] + "\n")
    else:
        dataout.append("FORMAMIDE_CONCENTRATION=0\n")
#    dataout.append("FORMAMIDE_CONC_UNITS=PERCENT\n") # Not supported
    
    if (g_omp_parameters["TMAC_CONCENTRATION"] != ""):
        dataout.append("TMAC_CONCENTRATION="+ g_omp_parameters["TMAC_CONCENTRATION"] + "\n")
    else:
        dataout.append("TMAC_CONCENTRATION=0\n")
    
    if (g_omp_parameters["BETAINE_CONCENTRATION"] != ""):
        dataout.append("BETAINE_CONCENTRATION=" + g_omp_parameters["BETAINE_CONCENTRATION"] + "\n")
    else:
        dataout.append("BETAINE_CONCENTRATION=0\n")
        
    if (g_omp_parameters["PH"] != ""):
        dataout.append("PH=" + g_omp_parameters["PH"] + "\n")
    else:
        dataout.append("PH=7\n")
    
    if (g_omp_parameters["POLYMER_SALT"] != ""):
        dataout.append("POLYMER_SALT=" + g_omp_parameters["POLYMER_SALT"] + "\n")
    else:
        dataout.append("POLYMER_SALT=TRUE\n")
    
    dataout.append("MICROCHIP_CORRECTION=No Corrections\n")
    dataout.append("SURFACE_SLOPE_DELG=1\n")
    dataout.append("SURFACE_INTERCEPT_DELG=0\n")
    dataout.append("SURFACE_SLOPE_DELH=1\n")
    dataout.append("SURFACE_INTERCEPT_DELH=0\n")
    
    dataout.append("\n[defaults]\n")
    dataout.append("STRAND_DEFAULT=Single\n")
    
    dataout.append("\n[sequences]\n")
    
    dataout.append("\nEXCLUDE_TGT_TGT_DUPLEX=TRUE\n")
    dataout.append("ILMAX_MONOMER=10\n")
    dataout.append("ILMAX_HOMODIMER=10\n")
    dataout.append("ILMAX_HETERODIMER=10\n")
    dataout.append("BULGEMAX_MONOMER=7\n")
    dataout.append("BULGEMAX_HOMODIMER=7\n")
    dataout.append("BULGEMAX_HETERODIMER=7\n")
    dataout.append("DEFAULT_GENERATESTRUCTURE_MONOMER=True\n")
    dataout.append("DEFAULT_GENERATESTRUCTURE_HOMODIMER=False\n")
    dataout.append("DEFAULT_GENERATESTRUCTURE_HETERODIMER=False\n")
    dataout.append("SUBOPTIMAL_ENABLE=TRUE\n")
    dataout.append("MAX_STRUCTURES_MONOMER=10\n")
    dataout.append("MAX_STRUCTURES_HOMODIMER=10\n")
    dataout.append("MAX_STRUCTURES_HETERODIMER=10\n")
#    dataout.append("SUBOPTIMALENERGY_WINDOW_ENABLE=TRUE\n") # Not supported
    dataout.append("SUBOPTIMALENERGY_WINDOW_MONOMER=1\n")
    dataout.append("SUBOPTIMALENERGY_WINDOW_HOMODIMER=1\n")
    dataout.append("SUBOPTIMALENERGY_WINDOW_HETERODIMER=1\n")
##    dataout.append("SUBOPTIMALWINDOW_WINDOW_ENABLE=TRUE\n") # Not supported
#    dataout.append("SUBOPTIMALWINDOW_MONOMER=90\n")
#    dataout.append("SUBOPTIMALWINDOW_HOMODIMER=90\n")
#    dataout.append("SUBOPTIMALWINDOW_HETERODIMER=90\n")
##    dataout.append("SUBOPTIMALDISTANCE_WINDOW_ENABLE=TRUE\n") # Not supported
#    dataout.append("SUBOPTIMALDISTANCE_MONOMER=7\n")
#    dataout.append("SUBOPTIMALDISTANCE_HOMODIMER=7\n")
#    dataout.append("SUBOPTIMALDISTANCE_HETERODIMER=7\n")

    dataout.append("\nSEQUENCE_NAME=" + probe_id + "\n")
    dataout.append("OLIGO_ENABLED=TRUE\n")
    dataout.append("SEQUENCE=" + probe_seq + "\n")
    
    if (g_omp_parameters["PROBE_CONCENTRATION"] != ""):
        dataout.append("CONCENTRATION=" + g_omp_parameters["PROBE_CONCENTRATION"] + "\n")
    else:
        dataout.append("CONCENTRATION=1E-07\n")
    
    dataout.append("SEQUENCE_TYPE=DNA\n")
    dataout.append("STRAND=SINGLE\n")
    dataout.append("SEQUENCE_FUNCTION=PROBE\n")
    dataout.append("TAIL_FOLDING=TRUE\n")
    dataout.append("EXCLUDE_SPECIES=" + probe_id + "," + probe_id + "\n")
    
    dataout.append("\nSEQUENCE_NAME=" + target_id + "\n")
    dataout.append("OLIGO_ENABLED=TRUE\n")
    dataout.append("SEQUENCE=" + target_seq + "\n")
    
    if (g_omp_parameters["TARGET_CONCENTRATION"] != ""):
        dataout.append("CONCENTRATION=" + g_omp_parameters["TARGET_CONCENTRATION"] + "\n")
    else:
        dataout.append("CONCENTRATION=1E-07\n")
    
    dataout.append("SEQUENCE_TYPE=" + target_type + "\n")
    dataout.append("STRAND=SINGLE\n")
    dataout.append("SEQUENCE_FUNCTION=TARGET\n")
    dataout.append("EXCLUDE_SPECIES_TYPE=" + target_id + ",TARGET\n")
    dataout.append("TAIL_FOLDING=TRUE\n")
    
#    print "======================================================================"
#    print "======================================================================"
#    print string.join(dataout, "")
#    print "======================================================================"
#    print "======================================================================"

    f = open(file_path, 'w')
    f.writelines(dataout)
    f.close()
    
def addLogEntry(log_line_list):
    log_file_path = g_omp_parameters["OMP_LOG_FILE_PATH"]
    log_file_found = os.path.exists(log_file_path)
    
    start_time = log_line_list[2]
    end_time = log_line_list[3]
    
    start_timestamp = start_time.strftime("%m/%d/%Y %H:%M:%S")
    end_timestamp = end_time.strftime("%m/%d/%Y %H:%M:%S")
    
    log_line_list[2] = str(start_timestamp)
    log_line_list[3] = str(end_timestamp)
    
    log_line = string.join(log_line_list, ",") + "\n"
    
    f = open(log_file_path, 'a')
    
    if (not log_file_found):
        header_list = ["Probe/Target Comparison File Path", "Number of Comparisons", 
                       "Start Time", "End Time", "Time to Completion (seconds)"]
        f.write(string.join(header_list, ",") + "\n")
    
    f.write(log_line)
    f.close()
    
if __name__ == "__main__":
    print "Running..."
    main()
    print "Done!!!"
