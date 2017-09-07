# particle_filter
A simple particle filter in 2D to track a vehicles position and heading based on GPS data and a known map.

The basic idea here to
* initialize the vehicle position and heading based on initial GPS data (with some sensor noise)
* predict the vehicle position and heading in the next time instant (50Hz) based on a simple bicycle model with non-zero yaw rate
* use sensor measurements onboard the vehicle to calculate likelihood weights for all the particles
* choose the most likely particles (resampling) and repeat 

Using this approach, over a few time steps, the particle filter weeds out unlikely positions for the vehicle and converges to the most likely position and heading for the moving vehicle

* 100 particles were chosen for this project - they satisfy the localization criteria well.
* A brute force nearest neighbor approach is used to associate map landmark data with the landmark measurements
