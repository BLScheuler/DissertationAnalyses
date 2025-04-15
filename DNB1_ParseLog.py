#!/usr/bin/python3
import os, re, gzip
DEBUG = False#True#

p_dict = {'obsidian': 1, 'redwood': 2, 'solstice': 3, 'starlight': 4, 'tundra': 5}
base_dir = '/Users/zjr109/git/Dissertation/RawData/DNB/'
#base_dir = '/Users/joehoupt/Documents/GitHub/Dissertation/RawData/DNB/'

islogfile = re.compile(".*\.log\.gz$")
participant_name = re.compile("(^\w+)_diss")

block_experimental = re.compile("^Instructions8Text.*text.*=")

block_remember_start = re.compile("^RightRemember.*autoDraw.*true")
block_remember_end = re.compile("^RightRemember.*autoDraw.*null")

block_break_start = re.compile("^MinBreakText.*autoDraw.*true")
block_break_end = re.compile("^MinBreak2Text.*autoDraw.*null")

block_AO = re.compile("^LeftClickAO.*numClicks.*null")
block_VO = re.compile("^LeftClickVO.*numClicks.*null")
block_EXP = re.compile("^LeftClickEXP.*numClicks.*null")
new_trial = re.compile("LeftClick.*numClicks.*null")
end_trial = re.compile("text.*autoDraw.*true")

highlightdoor = re.compile("HighlightDoor.*pos.*\[(\-?\d\.?\d?),")
sound = re.compile("Sound2.*secs.*2\.5")

mouse_down = re.compile("Mouse.*\s(\d)\sbutton\sdown")
mouse_up = re.compile("Mouse.*button\sup")
mouse_right = re.compile("RightClick.*numClicks.*=\s(\d+)$")
mouse_left = re.compile("LeftClick.*numClicks.*=\s(\d+)$")

dnb = open("dnb.csv", "w")
dnb.write("subject,week,block_no,trial_no,visual,auditory,choice,rt")
if DEBUG: 
  dnb.write(",condition\n")
else : 
  dnb.write("\n")

# Look for practGo3Pic for last screen before practice go (then 5)
# Look for GO_test_inst for last screen before test go (then 50)
# Look for PractStop3Pic for last screen before practice stop (then 5 x 2)
# Look for breakImage for last screen before test block

## Check abnormally long trial Meadow 46; extra trial in some blocks

for wk in range(3) : 
  files = os.listdir(base_dir + 'Week' + str(wk+1) + '/')
  for fname in files: 
    if islogfile.search(fname) :
      #fname = "Cascade_dissOSARI_2024-12-05_10h35.23.204.log"
      pname = participant_name.findall(fname)[0].lower()
      pnumber = p_dict[pname]
      press_onset = 0
      trial_state_last = 'null'
      trial_state = 'null'
      mouse_state = 'null'
      block_state_last = 'null'
      block_state = 'null'
      block_type_last = 'null'
      block_type = 'null'
      fill_state = 'null'
      fill_size = 0
      fill_size_last = 0
      fill_onset = 0
      fill_center = 0
      target_center = 0
      time_stamp_last = 0
      trial_n = 0
      block_n = 1
      state_stack = ''
      ready_experimental = -1
      trial_onset = 0
      door_cue = None
      sound_cue = None

      fname_full = base_dir + "Week" + str(wk + 1) + '/' + fname
      print(fname_full)
      with gzip.open(fname_full, 'rb') as infile:
        mouse_down_time = 0
        time_stamp = 0
        for line in infile:
          line = line.decode().strip()
          line = line.split('\t')
          if len(line) > 1 :
            time_stamp = round(float(line[0]),4)
  
            ## Read mouse state ##
            if line[1] == 'DATA':
              if mouse_down.search(line[2]): 
                press_onset = float(line[0])
                mouse_state = 'down'
                if trial_onset > 0 : 
                  rt = time_stamp - trial_onset
              elif mouse_up.search(line[2]): 
                press_duration = float(line[0]) - press_onset
                mouse_state = 'up'
  
            elif line[1] == 'EXP':
              ## Update Block State ##
              if block_experimental.search(line[2]):
                ready_experimental += 1
              elif block_remember_start.search(line[2]):
                block_state = 'remember'
              elif block_remember_end.search(line[2]):
                if ready_experimental > 0 :
                  block_state = "experimental"
              if block_AO.search(line[2]): 
                block_type = "ao"
              elif block_VO.search(line[2]): 
                block_type = "vo"
              elif block_EXP.search(line[2]): 
                block_type = "av"
              elif block_break_end.search(line[2]): 
                block_state = "break"

              if block_type != block_type_last : 
                  block_n = 1
   
                 ## 27 trials; look for remind or break..
              if block_state != block_state_last:
                if block_state == 'break': 
                  block_state = 'experimental'
                  block_n = block_n + 1
                  trial_n = 0
              
              if block_state == 'experimental' :
                ## Update Trial State ##
                if new_trial.search(line[2]): 
                  trial_state = 'begin'
                elif end_trial.search(line[2]):
                  trial_state = 'end'

                ## Set door cue
                if highlightdoor.search(line[2]):
                  door_cue = highlightdoor.findall(line[2])[0]
                elif sound.search(line[2]): 
                  sound_cue = '1'
  

                if mouse_right.search(line[2]) : 
                  if mouse_right.findall(line[2])[0] != '0' : 
                    choice = "right"
                if mouse_left.search(line[2]) : 
                  if mouse_left.findall(line[2])[0] != '0' : 
                    choice = "left"

                ## If trial state changed ##
                if trial_state != trial_state_last: 
                  # for testing, track state sequence 
                  if state_stack != '': 
                    state_stack += ","
                  state_stack += trial_state
  
                  if trial_state == 'begin' :  
                    door_cue = None
                    sound_cue = None
                    trial_type = None
                    rt = -1
                    trial_onset = time_stamp
                    trial_n = trial_n + 1
  
                  elif trial_state == 'end' :  
                    #dnb.write( str(pnumber) )
                    dnb.write(pname)
                    dnb.write(',' + str(wk+1) )
                    dnb.write(',' + str(block_n))
                    dnb.write(',' + str(trial_n))
                    if door_cue: 
                      dnb.write(',' + door_cue)
                    else: 
                      dnb.write(',None')
                    if sound_cue : 
                      dnb.write(',' + sound_cue)
                    else: 
                      dnb.write(',None')
                    dnb.write(',' + choice)
                    dnb.write(',' + str(round(rt,4)))
                    if DEBUG:
                      dnb.write(',' + block_type)
                      dnb.write(',' + str(time_stamp))
                      dnb.write(',' + pname)
                      #print(state_stack)
                    dnb.write("\n")


                    state_stack = ''
                
            trial_state_last = trial_state
            block_state_last = block_state
            block_type_last = block_type
            time_stamp_last = time_stamp
dnb.close()
