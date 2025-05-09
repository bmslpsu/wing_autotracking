# video2sparse

## how to use
`fly_work_dir`: the directory of the video folder   
`fly_folder`: video folder name  
`num_vid`: numbers of videos need to convert  
`start_frame`: the first frame to convert  
`end_frame`: the last frame to convert  
**make sure to have 15 frames before start frame and 15 frames after end_frame**  
`actionStartFrame`: any frame you want to be set as time 0  
`sparse_order312`: 1: switch into video order 312, 0: keep current video order. when testing, switching the video order to have the first video to be the one from a mirror.  
 
1. When running, choose the videos and backgrounds.  
2. Select the tether points from body to the end of tether.
3. Wait until done and check the converted files in folder start with: `Sparse_Frames_`