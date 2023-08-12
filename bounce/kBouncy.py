import numpy as np
import matplotlib.pyplot as plt
import time

def surface(x, a):
    return a * x ** 2

def surface_slope(x, a):
    return 2 * a * x

def update_ball_position(ball_position, ball_velocity, a, dt, t):
    ball_velocity[:, 1] += gravity * dt
    ball_position += ball_velocity * dt

    surface_y = surface(ball_position[:, 0], a)
    surface_bounce_mask = (ball_position[:, 1] <= surface_y)
    slope = surface_slope(ball_position[surface_bounce_mask, 0], a)
    normal = np.column_stack((-slope, np.ones_like(slope))) / np.sqrt(slope ** 2 + 1).reshape(-1, 1)
    ball_velocity[surface_bounce_mask] -= 2 * (ball_velocity[surface_bounce_mask] * normal).sum(axis=1).reshape(-1, 1) * normal
    ball_position[surface_bounce_mask, 1] = surface_y[surface_bounce_mask]

    return ball_position, ball_velocity

n_balls = 5000
a = 1
dt = 0.005
total_time = 60
n_frames = 1
frame_interval = int(total_time / (dt * n_frames))
x_max = 1
x = np.linspace(-1, 1, 2500)
y_surface = surface(x, a)

initial_x_positions = np.linspace(-1, 1, n_balls) * x_max
ball_positions = np.column_stack([initial_x_positions, np.full(n_balls, 1, dtype=float)])
ball_velocities = np.column_stack([np.full(n_balls, 0, dtype=float), np.full(n_balls, 0, dtype=float)])
gravity = -9.81
colors = plt.cm.turbo((initial_x_positions + x_max) / (2 * x_max))

frame_count = 0
previous_time = time.time()

frame_time = total_time / n_frames

time_step = 0
for t in np.arange(0, total_time, dt):
    update_ball_position(ball_positions, ball_velocities, a, dt, t)
    
    if frame_count < n_frames and time_step % frame_interval == 0:
        fig, ax1 = plt.subplots(1, 1, figsize=(6, 6), dpi=300)
        fig.tight_layout(pad=4.0)

        ax1.plot(x, y_surface, 'b-', alpha=0.5)
        for i in range(n_balls):
            ax1.plot(ball_positions[i][0], ball_positions[i][1], 'o', color=colors[i], markersize=1)

        ax1.set_xlim(-1, 1)
        ax1.set_ylim(0, 1)

        current_time = time.time()
        frame_time = current_time - previous_time
        mean_x = np.mean(ball_positions[:, 0])
        mean_y = np.mean(ball_positions[:, 1])

        ax1.text(0.475, 0.155, f'Time: {t:.2f} seconds', fontsize=5)
        ax1.text(0.395, 0.111, f'g: {(gravity):.2f}, num balls: {n_balls}', fontsize=5)
        ax1.text(0.335, 0.065, f'Avg. position: ({mean_x:.4f}, {mean_y:.4f})', fontsize=5)
        ax1.text(0.25, 0.02, f'Frame generation time: {frame_time:.2f}', fontsize=5)

        frame_count += 1
        plt.savefig(f'aa_{frame_count:05d}.png')
        plt.close(fig)
        previous_time = current_time
    time_step += 1

