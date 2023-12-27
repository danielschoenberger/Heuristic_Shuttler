import os
import tkinter as tk
from pathlib import Path

from PIL import Image, ImageTk


class ImageWidget:
    def __init__(self, frame_folder):
        self.frame_folder = frame_folder
        self.frame_files = sorted(os.listdir(frame_folder))
        self.num_frames = len(self.frame_files)

        self.root = tk.Tk()

        self.frame_index = tk.IntVar()
        self.frame_index.set(0)

        self.frame_slider = tk.Scale(
            self.root,
            from_=0,
            to=self.num_frames - 1,
            orient=tk.HORIZONTAL,
            length=200,
            variable=self.frame_index,
            command=self.update_image_from_slider,
        )

        self.buttons_frame = tk.Frame(self.root)
        self.prev_button = tk.Button(self.buttons_frame, text="Previous Frame", command=self.prev_frame)
        self.next_button = tk.Button(self.buttons_frame, text="Next Frame", command=self.next_frame)

        image_path = self.frame_folder / self.frame_files[0]
        self.pil_image = Image.open(image_path)
        self.tk_image = ImageTk.PhotoImage(self.pil_image)

        self.label = tk.Label(self.root, image=self.tk_image)
        self.label.pack()

        self.update_image()

    def prev_frame(self):
        self.frame_index.set((self.frame_index.get() - 1) % self.num_frames)
        self.update_image()

    def next_frame(self):
        self.frame_index.set((self.frame_index.get() + 1) % self.num_frames)
        self.update_image()

    def update_image_from_slider(self, _):
        self.update_image()

    def update_image(self):
        current_frame_index = self.frame_index.get()
        frame_path = self.frame_folder / self.frame_files[current_frame_index]
        self.pil_image = Image.open(frame_path)
        self.tk_image.paste(self.pil_image)
        self.label.config(image=self.tk_image)
        self.label.image = self.tk_image

    def display(self):
        self.frame_slider.pack()
        self.buttons_frame.pack()
        self.prev_button.pack(side=tk.LEFT)
        self.next_button.pack(side=tk.LEFT)
        self.root.mainloop()


# Get the absolute path to the script's directory
script_directory = Path(__file__).parent

# Replace 'your_folder_path' with the actual path to the folder containing your PNG frames
# Join the script's directory with the relative path to your folder
folder_path = script_directory / "run_20231222_130344"

# Create the ImageWidget instance
image_widget = ImageWidget(folder_path)

# Display the widget with slider, "Previous Frame" button, and "Next Frame" button below the slider
image_widget.display()
