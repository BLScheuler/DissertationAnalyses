#!/usr/bin/python3
import os, re, gzip
DEBUG = False

p_dict = {'cascade': 1, 'glacier': 2, 'harbor': 3, 'horizon': 4, 'meadow': 5}
#base_dir = '/Users/zjr109/git/Dissertation/RawData/OSARI/'
base_dir = '/Users/joehoupt/Documents/GitHub/Dissertation/RawData/OSARI/'

islogfile = re.compile(".*\.log\.gz$")
participant_name = re.compile("(^\w+)_diss")

block_pracgo = re.compile("Go_inst:.*image")
#block_break = re.compile("breakImage.*autodraw.*=\s(\w)$")
block_break = re.compile("breakImage.*autoDraw\s*=\s+(\w*)$")
block_go = re.compile("GO_test_inst.*autoDraw\s*=\s+(\w*)$")
block_endgo =re.compile("mixIntro.*autoDraw\s*=\s+(\w*)$") 
block_mixprac = re.compile("PractStop3Pic.*autoDraw\s*=\s+(\w*)$")

trial_feedback = re.compile("feedbackTxt")
#trial_feedback = re.compile("FeedbackTxt.*text\s+=\s*")
#new_trial = re.compile("targetTriangle_left_2.*null")
timer = re.compile("timerText.*=\s(\d*\.?\d*)\s*$")
new_trial = re.compile("fillingLine.*size.*,0]$")

mouse_down = re.compile("Mouse.*button\sdown")
mouse_up = re.compile("Mouse.*button\sup")
fill_line = re.compile("fillingLine.*size.*,(\d+\.?\d*)]")
fill_pos = re.compile("fillingLine.*pos.*,(\d+\.?\d*)]")
fill_color = re.compile("fillColor\s=\s(.*)$")
targ_pos = re.compile("targetTriangle_left.*pos.*,(\d+\.?\d*)]")

osari = open("osari.csv", "w")
#osari.write("subject,click_duration,fill_y,target_y,trial_type\n")
osari.write("subject,ss_presented,inhibited,ssd,rt,week")
if DEBUG: 
  osari.write(",trial_no,block_no,condition\n")
else : 
  osari.write("\n")

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
      fill_state = 'null'
      fill_size = 0
      fill_size_last = 0
      fill_onset = 0
      fill_center = 0
      target_center = 0
      time_stamp_last = 0
      trial_n = 0
      block_n = 0
      state_stack = ''

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
              elif mouse_up.search(line[2]): 
                press_duration = float(line[0]) - press_onset
                mouse_state = 'up'
  
            elif line[1] == 'EXP':
              ## Update Block State ##
              if block_go.search(line[2]):
                if block_go.findall(line[2])[0]=='true' : 
                  block_state = 'go_inst'
                elif block_go.findall(line[2])[0]=='null' and block_state_last == 'go_inst':
                  block_state = 'go'
              elif block_endgo.search(line[2]):
                if block_endgo.findall(line[2])[0]=='true' : 
                  block_state = 'null'
              elif block_mixprac.search(line[2]):
                if block_mixprac.findall(line[2])[0]=='true' : 
                  block_state = 'mixprac_inst'
                elif block_mixprac.findall(line[2])[0]=='null' and block_state_last == 'mixprac_inst': 
                  block_state = 'mix_prac'
              elif block_break.search(line[2]):
                if block_break.findall(line[2])[0]=='true' : 
                  block_state = 'break'
                elif block_break.findall(line[2])[0]=='null' and block_state_last == 'break':
                    block_state = 'mix'
   
              if block_state != block_state_last:
                if block_state == 'break' or block_state=='go_inst': 
                  block_n = block_n + 1
                  trial_n = 0
              
              if block_state == 'go' or block_state == 'mix' : 
                ## Update Trial State ##
                if new_trial.search(line[2]): 
                  trial_state = 'begin'
                elif fill_line.search(line[2]) : 
                  fill_size = round(float(fill_line.findall(line[2])[0]), 3)
                  if fill_size != fill_size_last : 
                    trial_state = 'filling'
                  elif trial_state_last == 'filling': 
                    trial_state = 'stopped'
                elif trial_feedback.search(line[2]): 
                  trial_state = 'feedback'
  
  
                ## If trial state changed ##
                if trial_state != trial_state_last: 
                  # for testing, track state sequence 
                  if state_stack != '': 
                    state_stack += ","
                  state_stack += trial_state
  
                  if trial_state == 'begin' :  
                    press_onset = 0
                    press_duration = 0
                    fill_onset = 0
                    fill_y = 0
                    trial_n = trial_n + 1
                    trial_type = "go" # Always start as go, set to stop if stopped filling
                    fill_state = 'empty'
  
                  elif trial_state == 'filling' :  
                    fill_onset = time_stamp_last
                    if press_onset == 0 :   # Deal with held down button after instructions
                      press_onset = time_stamp_last
  
                  elif trial_state == 'stopped' :  
                    if mouse_state == 'down' and fill_y < .14:
                      trial_type="stop"
                      ssd = time_stamp_last - fill_onset
  
                  elif trial_state == 'feedback' :  
                    #osari.write( str(3*(pnumber-1) + wk + 1))
                    osari.write( str(pnumber) )
                    if trial_type == "stop": 
                      osari.write(',1')
                      if press_duration ==1 or arrow_color == "#008000":
                        osari.write(',1,500')
                      else: 
                        osari.write(',0,500')
                    else :
                      osari.write(',0')
                      osari.write(',-999,-999')
                    osari.write(',' + str(round(1000*press_duration,0)))
                    osari.write(',' + str(wk+1))
                    if DEBUG:
                      osari.write(',' + str(trial_n))
                      osari.write(',' + str(block_n))
                      osari.write(',' + block_state)
                      osari.write(',' + str(time_stamp))
                      osari.write(',' + pname)
                      print(state_stack)
                    osari.write("\n")
  
                    state_stack = ''
                
                if timer.search(line[2]) :
                  if float(timer.findall(line[2])[0]) >= 1:
                    press_duration = 1
  
                if fill_pos.search(line[2]) : 
                  fill_center = round(float(fill_pos.findall(line[2])[0]), 3)
                  fill_y = fill_center + fill_size / 2
  
                if fill_color.search(line[2]) : 
                    arrow_color = fill_color.findall(line[2])[0]
  
            fill_size_last = fill_size
            trial_state_last = trial_state
            block_state_last = block_state
            time_stamp_last = time_stamp
osari.close()
