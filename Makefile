OBJS = Across_pts.o Geocoordinate.o across_point_atmosphere.o beta.o \
	calculate_intensity.o convert_coordinate.o \
	cross_point_at_altitude.o \
	get_arg.o get_observation_data.o intensity_integral.o \
	interpolated_m.o \
	main.o obtainfittingcoefficient.o obtainfittingpmc.o \
	optical_depth.o \
	pmc.o search_pmc_region.o select_latlon.o solar_direction.o \
	tangential_point.o this_position_is_near_pmc.o does_pmc_exist.o \
	output.o write_log.o

LIBS = -L/home/ando/lib -Wl,-rpath /home/ando/lib -lnlopt -lMie_scattering -lcomplex_bessel -lnrlmsise21 -lAndoLab_22
INCS = -I/home/ando/include
OPTS = -Wall -O3 -std=c++17 $(INCS)

.PHONY: all clean

all: main

main: $(OBJS)
	g++ -o $@ $(OBJS) $(LIBS)

%.o: %.cpp
	g++ -c $< $(OPTS)

intensity_integra.o optical_depth.o Across_pts.o calculate_intensity.o search_pmc_region.o: Across_pts.h

intensity_integral.o optical_depth.o get_arg.o main.o solar_direction.o: Date.h

select_latlon.o get_observation_data.o intensity_integral.o Geocoordinate.o across_point_atmosphere.o Across_pts.o calculate_intensity.o main.o solar_direction.o pmc.o: Fixed_parameter.h

intensity_integral.o optical_depth.o beta.o Geocoordinate.o main.o: Geocoordinate.h

optical_depth.o search_pmc_region.o: Region.h

optical_depth.o main.o solar_direction.o pmc.o: pmc.h

select_latlon.o get_observation_data.o obtainfittingcoefficient.o main.o: pmc_observation.h

tangential_point.o convert_coordinate.o intensity_integra.o optical_depth.o cross_point_at_altitude.o beta.o across_point_atmosphere.o obtainfittingcoefficient.o obtainfittingpmc.o calculate_intensity.o get_arg.o this_position_is_near_pmc.o main.o search_pmc_region.o: pmc_simulation.h

main.o write_log.o: write_log.h

