'''
Python module attempting to simulate the formation of a solar system from randomly dispersed particles
with randomly initial velocities under the influence of gravity.
'''
import collections
import random
import tqdm
import math
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import matplotlib.cm as cm
import traceback

import cython_force_calc
import cython_collision

Position = collections.namedtuple('Position', ['x', 'y', 'z'])
Velocity = collections.namedtuple('Velocity', ['x', 'y', 'z'])



class SolarSim():

    TIME_STEP_DAYS = 5
    TIME_STEP_SECONDS = TIME_STEP_DAYS * 24 * 60 * 60

    CAMERA_ROTATION_RATE = 0
    SHOW_ANIMATION = False

    MAX_DISTANCE_OF_IMPACT = 10 ** 11.1  # 11.0
    GRAVITATIONAL_CONSTANT = 6.672 * 10 ** (-11)  # N*m^2*kg^-2
    TOTAL_MASS = 2 * (1.9895 * 10 ** 30)  # 2X mass of the solar system (a lot of the initial mass drifts off/gets expelled.  Starting 2x diameter, so 2x mass seems reasonable)

    MAX_COLOR_VALUE = 10 ** 2.6
    MASS_TO_SIZE_RATIO = 1 / 10 ** 17.6

    # 10**12.8 meters is roughly diameter of our solar system (Pluto's orbit)
    MAX_X_INITIAL_POS = 2 * (10 ** 12.8)  # meters  --> starting 2x solar system diameter
    MAX_Y_INITIAL_POS = 2 * (10 ** 12.8)  # meters --> starting 2x solar system diameter
    MAX_Z_INITIAL_POS = 10 ** 12.5  # meters  # was 10**12.3
    BOX_DIAGONAL_DIST = (MAX_X_INITIAL_POS**2 + MAX_Y_INITIAL_POS**2 + MAX_Z_INITIAL_POS**2)**0.5
    REMOVAL_DISTANCE = BOX_DIAGONAL_DIST * 0.9  # 1.0
    CENTER_X = MAX_X_INITIAL_POS / 2
    CENTER_Y = MAX_Y_INITIAL_POS / 2
    CENTER_Z = MAX_Z_INITIAL_POS / 2
    MAX_INITIAL_ABS_VEL = 10 ** 3.35  # meters/second
    MAX_ROTATE_VEL = 10 ** 4.0


    def __init__(self, TOTAL_TIME=10000, INITIAL_SPACEBODIES=500, FRAME_SAMPLE_RATE=10, ANIMATION_INTERVAL=30, ANIMATION_FILENAME_LETTER=None):
        self.TOTAL_TIME = TOTAL_TIME
        self.NUMBER_OF_INITIAL_SPACEBODIES = INITIAL_SPACEBODIES
        self.AVG_INITIAL_MASS = self.TOTAL_MASS / self.NUMBER_OF_INITIAL_SPACEBODIES
        self.FRAME_SAMPLE_RATE = FRAME_SAMPLE_RATE
        self.ANIMATION_FILENAME_LETTER = ANIMATION_FILENAME_LETTER
        self.ANIMATION_INTERVAL = ANIMATION_INTERVAL
        if ANIMATION_FILENAME_LETTER is None:
            self.ANIMATION_FILENAME = f'SimulationVideos\\SolarSim - TOTAL_TIME {self.TOTAL_TIME}  NUM_SB {self.NUMBER_OF_INITIAL_SPACEBODIES}'
        else:
            self.ANIMATION_FILENAME = f'SimulationVideos\\SolarSim - {self.ANIMATION_FILENAME_LETTER} - TOTAL_TIME {self.TOTAL_TIME}  NUM_SB {self.NUMBER_OF_INITIAL_SPACEBODIES}'

        self.solar_system_data = self.simulate_solar_system()
        self.solar_system_data = {period: spacebodies for period, spacebodies in self.solar_system_data.items() if period % self.FRAME_SAMPLE_RATE == 0}

        self.solar_system_animation(ELEVATION_ANGLE=60)
        self.solar_system_animation(ELEVATION_ANGLE=45)
        self.solar_system_animation(ELEVATION_ANGLE=0)


    def solar_system_animation(self, ELEVATION_ANGLE=45):
        fig = plt.figure(figsize=(12, 8))
        ax = p3.Axes3D(fig)
        title = ax.set_title('Solar System Simulation   -   Period: 0')
        plt.setp(title, color=(0.8, 0.8, 1.0, 1.0))
        ax.patch.set_facecolor((0.05, 0.05, 0.05, 1.0))
        ax.set_axis_off()
        ax.view_init(elev=ELEVATION_ANGLE)

        ax.set_xlim3d([0, self.MAX_X_INITIAL_POS])
        ax.set_xlabel('X')
        ax.set_ylim3d([0, self.MAX_Y_INITIAL_POS])
        ax.set_ylabel('Y')
        ax.set_zlim3d([0, self.MAX_Z_INITIAL_POS])
        ax.set_zlabel('Z')
        ax.set_zlim([self.CENTER_Z - self.MAX_X_INITIAL_POS / 2, self.CENTER_Z + self.MAX_X_INITIAL_POS / 2])

        initial_x = [sb.position[0] for sb in self.solar_system_data[0]]
        initial_y = [sb.position[1] for sb in self.solar_system_data[0]]
        initial_z = [sb.position[2] for sb in self.solar_system_data[0]]
        initial_size = [area_from_vol(sb.mass) * self.MASS_TO_SIZE_RATIO for sb in self.solar_system_data[0]]

        cmap = cm.get_cmap('plasma')
        normalize = matplotlib.colors.Normalize(vmin=0, vmax=self.MAX_COLOR_VALUE)
        initial_colors = [cmap(normalize(size)) for size in initial_size]
        scatter = ax.scatter(initial_x, initial_y, initial_z, s=initial_size, c=initial_colors, vmin=0, vmax=self.MAX_COLOR_VALUE, alpha=1.0)

        pbar = tqdm.tqdm(total=self.TOTAL_TIME)

        def update(frame_number):
            x, y, z = list(zip(*[(sb.position[0], sb.position[1], sb.position[2]) for sb in self.solar_system_data[frame_number]]))  # fancy instead of 3 listcomps
            size = [area_from_vol(sb.mass) * self.MASS_TO_SIZE_RATIO for sb in self.solar_system_data[frame_number]]
            colors = [cmap(normalize(s)) for s in size]

            scatter._offsets3d = (x, y, z)
            scatter.set_sizes(size)
            scatter._facecolor3d = colors  # HAVE TO OVERWRITE THE COLOR ATTRIBUTES AS A MATPLOTLIB BUG IS MESSING THEM UP!!!
            scatter._edgecolor3d = colors  # HAVE TO OVERWRITE THE COLOR ATTRIBUTES AS A MATPLOTLIB BUG IS MESSING THEM UP!!!
            title.set_text(f'Solar System Simulation   -   Period:  {frame_number}')
            if self.CAMERA_ROTATION_RATE != 0:
                ax.azim += self.CAMERA_ROTATION_RATE
            if frame_number % 100 == 0:
                pbar.update(100)


        ss_animation = animation.FuncAnimation(fig, update, frames=[period for period in range(self.TOTAL_TIME) if period % self.FRAME_SAMPLE_RATE == 0], interval=self.ANIMATION_INTERVAL)
        # NEED FFMPEG INSTALLED FOR ANIMATION SAVING/CONVERSION TO .MP4!!!
        try:
            ss_animation.save(self.ANIMATION_FILENAME + f'  ELEV_ANGLE {ELEVATION_ANGLE}' + '.mp4')
        except:
            traceback.print_exc()
            print('ERROR!  Unable to save animation.')
            print('Please check that ffmpeg is installed and added to your PATH!')

        # See:  https://matplotlib.org/api/_as_gen/matplotlib.pyplot.close.html
        plt.close(fig)  # clears memory of fig object; otherwise memory can overflow
        if self.SHOW_ANIMATION:
            plt.show()


    def simulate_solar_system(self) -> dict:
        spacebodies = self.create_initial_environment()
        solar_system_data = {0: [SpaceBody(sb.id, sb.mass, sb.position, sb.velocity) for sb in spacebodies]}
        pbar = tqdm.tqdm(total=self.TOTAL_TIME)
        for period in range(1, self.TOTAL_TIME):
            spacebodies_data = [(sb.id, sb.mass, sb.position.x, sb.position.y, sb.position.z) for sb in spacebodies]
            for sb in spacebodies:
                sb.set_new_velocity(spacebodies_data)
                sb.set_new_position()

            spacebodies_copy = [SpaceBody(sb.id, sb.mass, sb.position, sb.velocity) for sb in spacebodies]
            spacebodies = cython_collision.calc_collisions(spacebodies,
                                                                          spacebodies_copy,
                                                                          self.TOTAL_MASS,
                                                                          self.MAX_DISTANCE_OF_IMPACT)

            if period % self.FRAME_SAMPLE_RATE == 0:
                solar_system_data[period] = [SpaceBody(sb.id, sb.mass, sb.position, sb.velocity) for sb in spacebodies]  # MUUUUCH faster with manual copy instead of copy.deepcopy!!!
                pbar.update(self.FRAME_SAMPLE_RATE)

        pbar.update(self.TOTAL_TIME % self.FRAME_SAMPLE_RATE)
        return solar_system_data


    def create_initial_environment(self) -> list:
        spacebodies = []

        def angle(x1: float, x2: float, y1: float, y2: float) -> float:
            '''
            Calcs angle (0-360) between 2 points (relative to the first point)
            with both points' x and y coordinates.
            '''
            dx, dy = x2 - x1, y2 - y1
            h = (dy ** 2 + dx ** 2) ** 0.5
            theta = math.degrees(math.acos(dx / h))
            if dy < 0:
                theta = math.degrees(math.acos(-dx / h))
                theta += 180
            return theta

        for id in range(self.NUMBER_OF_INITIAL_SPACEBODIES):
            # Including a lower bound on mass > 0 so we don't have unseeably small points that have little effect other than to slow down the simulation.
            initial_mass = random.uniform(self.AVG_INITIAL_MASS * 0.25, self.AVG_INITIAL_MASS * 1.75)

            initial_x_pos = random.gauss(self.MAX_X_INITIAL_POS / 2, self.MAX_X_INITIAL_POS * 0.25)  # gaussian dist
            initial_y_pos = random.gauss(self.MAX_Y_INITIAL_POS / 2, self.MAX_Y_INITIAL_POS * 0.25)  # gaussian dist
            initial_z_pos = random.gauss(self.MAX_Z_INITIAL_POS / 2, self.MAX_Z_INITIAL_POS * 0.25)  # gaussian dist

            initial_x_vel = random.uniform(-self.MAX_INITIAL_ABS_VEL, self.MAX_INITIAL_ABS_VEL)
            initial_y_vel = random.uniform(-self.MAX_INITIAL_ABS_VEL, self.MAX_INITIAL_ABS_VEL)
            initial_z_vel = random.uniform(-self.MAX_INITIAL_ABS_VEL, self.MAX_INITIAL_ABS_VEL)

            new_sb = SpaceBody(str(id),
                                         initial_mass,
                                         Position(initial_x_pos, initial_y_pos, initial_z_pos),
                                         Velocity(initial_x_vel, initial_y_vel, initial_z_vel))

            # ADDING ANGULAR VELOCITY
            xy_angle_to_center = angle(self.CENTER_X, initial_x_pos, self.CENTER_Y, initial_y_pos)
            xy_angle_to_rotate = xy_angle_to_center - 90

            dist_weight = new_sb.distance_from_center() / self.BOX_DIAGONAL_DIST

            x_vel_change = math.cos(math.radians(xy_angle_to_rotate)) * self.MAX_ROTATE_VEL * dist_weight * random.uniform(0.6, 1)
            y_vel_change = math.sin(math.radians(xy_angle_to_rotate)) * self.MAX_ROTATE_VEL * dist_weight * random.uniform(0.6, 1)
            new_velocity = Velocity(new_sb.velocity.x + x_vel_change,
                                                 new_sb.velocity.y + y_vel_change,
                                                 new_sb.velocity.z)
            new_sb.velocity = new_velocity

            spacebodies.append(new_sb)
        return spacebodies



class SpaceBody():

    def __init__(self, id: str, mass: float, initial_position, initial_velocity_vector) -> None:
        self.id = id
        self.mass = mass
        self.t_mass = SolarSim.TIME_STEP_SECONDS / mass
        self.position = initial_position
        self.velocity = initial_velocity_vector


    def set_new_velocity(self, spacebodies_data):
        # g_force = cython_force_calc.calc_composite_force_vector(SolarSim.GRAVITATIONAL_CONSTANT, self, all_spacebodies_list)
        g_force_x, g_force_y, g_force_z = cython_force_calc.calc_composite_force_vector(SolarSim.GRAVITATIONAL_CONSTANT, self, spacebodies_data)

        # V2 = V1 + a*t    -->   a = Force/Mass
        self_velocity, t_mass = self.velocity, self.t_mass  # local variables for speed

        self.velocity = Velocity(self_velocity[0] + g_force_x * t_mass,
                                           self_velocity[1] + g_force_y * t_mass,
                                           self_velocity[2] + g_force_z * t_mass)


    def set_new_position(self):
        self_position, self_velocity = self.position, self.velocity  # local variables for speed
        time_step_seconds = SolarSim.TIME_STEP_SECONDS  # local variable for speed
        self.position = Position(self_position[0] + self_velocity[0] * time_step_seconds,
                                           self_position[1] + self_velocity[1] * time_step_seconds,
                                           self_position[2] + self_velocity[2] * time_step_seconds)


    def distance_from_center(self) -> float:
        return ((self.position[0] - SolarSim.CENTER_X) ** 2 + (self.position[1] - SolarSim.CENTER_Y) ** 2 + (self.position[2] - SolarSim.CENTER_Z) ** 2) ** 0.5


    def __add__(self, other):
        new_id = f'{self.id}+{other.id}'
        new_mass = self.mass + other.mass
        self_pct, other_pct = self.mass / new_mass, other.mass / new_mass
        new_position = Position(self.position.x * (self_pct) + other.position.x * (other_pct),
                                             self.position.y * (self_pct) + other.position.y * (other_pct),
                                             self.position.z * (self_pct) + other.position.z * (other_pct))
        new_velocity = Velocity(self.velocity.x * (self_pct) + other.velocity.x * (other_pct),
                                             self.velocity.y * (self_pct) + other.velocity.y * (other_pct),
                                             self.velocity.z * (self_pct) + other.velocity.z * (other_pct))
        return SpaceBody(new_id, new_mass, new_position, new_velocity)


    def __repr__(self):
        return f'SpaceBody Object {self.id if self.id is not None else "NO ID"}'



def area_from_vol(volume: float) -> float:
    '''
    Takes a given Volume (saying mass ~= volume for SolarSim) and
    returns the equivalent area for the same radius.
    This is used b/c matplotlib plots points with sizes as AREA.
    Converting SpaceBody.mass from a "volume" to an area makes the points visually scale correctly.

    Math:

    Area = pi*radius**2
    Volume = (4/3)*pi*r**3

    Radius_from_Volume
    r = (Volume / ((4/3)*pi))**(1/3)

    Area from Volume
    Area = pi * (r) **2
    Area = pi * (Volume / ((4/3)*pi)) ** (2/3)
    '''
    return math.pi * (volume / 4.1888) ** (0.6667)



if __name__ == '__main__':
    x = SolarSim(TOTAL_TIME=10000, INITIAL_SPACEBODIES=150, ANIMATION_INTERVAL=50, FRAME_SAMPLE_RATE=15)
