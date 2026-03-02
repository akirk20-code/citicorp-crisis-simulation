"""
Citicorp Center — OpenCV Point Picking for Structural Dimension Extraction
==========================================================================
Click on structural points in a construction photo or diagram to extract
real dimensions, angles, and proportions.

Usage:
    python extract_dimensions.py [image_path]

Controls:
    Left click   — Pick a point
    Right click  — Undo last point
    'r'          — Set reference (pick 2 points, enter known distance)
    'a'          — Compute angle between last 3 points
    'd'          — Compute distance between last 2 points
    'm'          — Measure all consecutive pairs
    'e'          — Export results to CSV
    'c'          — Clear all points
    'z'          — Zoom toggle (click to set zoom center)
    'q' / ESC    — Quit

Author: Generated with Claude Code (OR 750, GMU)
"""

import cv2
import numpy as np
import sys
import csv
import os
from datetime import datetime


class PointPicker:
    def __init__(self, image_path):
        self.img_orig = cv2.imread(image_path)
        if self.img_orig is None:
            raise FileNotFoundError(f"Cannot load image: {image_path}")

        self.image_path = image_path
        self.img = self.img_orig.copy()
        self.points = []           # List of (x, y) pixel coordinates
        self.labels = []           # Optional label for each point
        self.ref_scale = None      # pixels per ft (set with 'r')
        self.ref_dist_ft = None    # reference distance in ft
        self.zoom_active = False
        self.zoom_center = None
        self.zoom_factor = 3.0
        self.window_name = "Citicorp Center — Point Picker"

        # Results
        self.measurements = []  # List of dicts: {type, points, value, unit}

        cv2.namedWindow(self.window_name, cv2.WINDOW_NORMAL)
        cv2.resizeWindow(self.window_name, 1200, 800)
        cv2.setMouseCallback(self.window_name, self._mouse_cb)

    def _mouse_cb(self, event, x, y, flags, param):
        if event == cv2.EVENT_LBUTTONDOWN:
            # Convert from display coords to image coords if zoomed
            ix, iy = self._display_to_image(x, y)
            self.points.append((ix, iy))

            # Ask for label
            label = f"P{len(self.points)}"
            self.labels.append(label)
            print(f"  Point {label}: ({ix}, {iy}) px", end="")
            if self.ref_scale is not None:
                print(f"  [{ix/self.ref_scale:.1f}, {iy/self.ref_scale:.1f}] ft", end="")
            print()
            self._redraw()

        elif event == cv2.EVENT_RBUTTONDOWN:
            if self.points:
                removed = self.points.pop()
                self.labels.pop()
                print(f"  Removed point at ({removed[0]}, {removed[1]})")
                self._redraw()

    def _display_to_image(self, dx, dy):
        """Convert display coordinates to image coordinates (accounts for zoom)."""
        if self.zoom_active and self.zoom_center is not None:
            # Get display window size
            h, w = self.img_orig.shape[:2]
            zh = int(h / self.zoom_factor)
            zw = int(w / self.zoom_factor)
            cx, cy = self.zoom_center
            x0 = max(0, cx - zw // 2)
            y0 = max(0, cy - zh // 2)

            # Scale display coords back to image coords
            win_w, win_h = cv2.getWindowImageRect(self.window_name)[2:]
            if win_w > 0 and win_h > 0:
                ix = x0 + int(dx * zw / win_w)
                iy = y0 + int(dy * zh / win_h)
                return ix, iy
        return dx, dy

    def _redraw(self):
        """Redraw image with all picked points and connections."""
        self.img = self.img_orig.copy()

        # Draw points
        for i, (px, py) in enumerate(self.points):
            color = (0, 255, 0) if i < len(self.points) - 1 else (0, 0, 255)
            cv2.circle(self.img, (px, py), 5, color, -1)
            cv2.circle(self.img, (px, py), 7, (255, 255, 255), 1)
            # Label
            label = self.labels[i] if i < len(self.labels) else f"P{i+1}"
            cv2.putText(self.img, label, (px + 10, py - 10),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 0), 1)

        # Draw lines between consecutive points
        for i in range(len(self.points) - 1):
            p1 = self.points[i]
            p2 = self.points[i + 1]
            cv2.line(self.img, p1, p2, (0, 200, 200), 1, cv2.LINE_AA)

        # Status bar
        h, w = self.img.shape[:2]
        bar_h = 40
        cv2.rectangle(self.img, (0, h - bar_h), (w, h), (40, 40, 40), -1)
        status = f"Points: {len(self.points)}  |  "
        if self.ref_scale is not None:
            status += f"Scale: {self.ref_scale:.2f} px/ft  |  "
        else:
            status += "Scale: not set (press 'r')  |  "
        status += "LClick: pick | RClick: undo | r: ref | d: dist | a: angle | e: export | q: quit"
        cv2.putText(self.img, status, (10, h - 12),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.4, (200, 200, 200), 1)

        # Apply zoom if active
        display = self.img
        if self.zoom_active and self.zoom_center is not None:
            zh = int(h / self.zoom_factor)
            zw = int(w / self.zoom_factor)
            cx, cy = self.zoom_center
            x0 = max(0, min(cx - zw // 2, w - zw))
            y0 = max(0, min(cy - zh // 2, h - zh))
            display = self.img[y0:y0+zh, x0:x0+zw]

        cv2.imshow(self.window_name, display)

    def set_reference(self):
        """Use the last 2 points to set the pixel-to-ft scale."""
        if len(self.points) < 2:
            print("  Need at least 2 points to set reference.")
            return

        p1 = self.points[-2]
        p2 = self.points[-1]
        px_dist = np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

        print(f"  Reference points: {p1} -> {p2}")
        print(f"  Pixel distance: {px_dist:.1f} px")
        try:
            ft_str = input("  Enter known distance in ft: ").strip()
            ft_dist = float(ft_str)
            self.ref_scale = px_dist / ft_dist
            self.ref_dist_ft = ft_dist
            self.measurements.append({
                'type': 'reference',
                'points': [p1, p2],
                'px': px_dist,
                'ft': ft_dist,
                'scale': self.ref_scale
            })
            print(f"  Scale set: {self.ref_scale:.2f} px/ft ({ft_dist} ft = {px_dist:.1f} px)")
        except ValueError:
            print("  Invalid input. Scale not set.")

    def measure_distance(self):
        """Measure distance between the last 2 points."""
        if len(self.points) < 2:
            print("  Need at least 2 points.")
            return

        p1 = self.points[-2]
        p2 = self.points[-1]
        px_dist = np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

        result = {
            'type': 'distance',
            'points': [p1, p2],
            'labels': [self.labels[-2], self.labels[-1]],
            'px': px_dist
        }

        print(f"  Distance {self.labels[-2]}->{self.labels[-1]}: {px_dist:.1f} px", end="")
        if self.ref_scale is not None:
            ft_dist = px_dist / self.ref_scale
            result['ft'] = ft_dist
            print(f" = {ft_dist:.2f} ft", end="")
        print()

        # Horizontal and vertical components
        dx = abs(p2[0] - p1[0])
        dy = abs(p2[1] - p1[1])
        print(f"    dx={dx:.1f} px, dy={dy:.1f} px", end="")
        if self.ref_scale is not None:
            print(f"  ({dx/self.ref_scale:.2f} ft, {dy/self.ref_scale:.2f} ft)", end="")
        print()

        self.measurements.append(result)

    def measure_angle(self):
        """Measure angle at the middle of the last 3 points (vertex = middle point)."""
        if len(self.points) < 3:
            print("  Need at least 3 points (vertex = middle point).")
            return

        p1 = np.array(self.points[-3])
        p2 = np.array(self.points[-2])  # vertex
        p3 = np.array(self.points[-1])

        v1 = p1 - p2
        v2 = p3 - p2

        cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
        angle_deg = np.degrees(np.arccos(np.clip(cos_a, -1, 1)))

        # Also compute angle from horizontal for each leg
        angle_h1 = np.degrees(np.arctan2(abs(v1[1]), abs(v1[0])))
        angle_h2 = np.degrees(np.arctan2(abs(v2[1]), abs(v2[0])))

        result = {
            'type': 'angle',
            'points': [tuple(p1), tuple(p2), tuple(p3)],
            'labels': [self.labels[-3], self.labels[-2], self.labels[-1]],
            'angle_deg': angle_deg,
            'leg1_from_horiz': angle_h1,
            'leg2_from_horiz': angle_h2
        }

        print(f"  Angle at {self.labels[-2]}: {angle_deg:.1f} deg")
        print(f"    Leg {self.labels[-3]}-{self.labels[-2]}: {angle_h1:.1f} deg from horizontal")
        print(f"    Leg {self.labels[-2]}-{self.labels[-1]}: {angle_h2:.1f} deg from horizontal")

        self.measurements.append(result)

        # Draw angle arc on image
        self._draw_angle_arc(tuple(p2), tuple(p1), tuple(p3), angle_deg)

    def _draw_angle_arc(self, vertex, p1, p3, angle_deg):
        """Draw an angle arc on the image."""
        radius = 30
        v1 = np.array(p1) - np.array(vertex)
        v2 = np.array(p3) - np.array(vertex)
        a1 = np.degrees(np.arctan2(-v1[1], v1[0]))  # note: y is flipped in image coords
        a2 = np.degrees(np.arctan2(-v2[1], v2[0]))
        cv2.ellipse(self.img, vertex, (radius, radius), 0, -a1, -a2, (0, 255, 255), 2)
        # Label
        mid_angle = (a1 + a2) / 2
        lx = int(vertex[0] + (radius + 15) * np.cos(np.radians(mid_angle)))
        ly = int(vertex[1] - (radius + 15) * np.sin(np.radians(mid_angle)))
        cv2.putText(self.img, f"{angle_deg:.1f}", (lx, ly),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 255), 1)
        cv2.imshow(self.window_name, self.img)

    def measure_all_pairs(self):
        """Measure distances between all consecutive point pairs."""
        if len(self.points) < 2:
            print("  Need at least 2 points.")
            return

        print(f"\n  === All consecutive distances ({len(self.points)-1} segments) ===")
        for i in range(len(self.points) - 1):
            p1 = self.points[i]
            p2 = self.points[i + 1]
            px_dist = np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)
            dx = abs(p2[0] - p1[0])
            dy = abs(p2[1] - p1[1])
            angle_h = np.degrees(np.arctan2(dy, dx))

            line = f"  {self.labels[i]:>4s} -> {self.labels[i+1]:<4s}: {px_dist:7.1f} px"
            if self.ref_scale is not None:
                ft_dist = px_dist / self.ref_scale
                dx_ft = dx / self.ref_scale
                dy_ft = dy / self.ref_scale
                line += f" = {ft_dist:7.2f} ft  (dx={dx_ft:.1f}, dy={dy_ft:.1f} ft)"
            line += f"  angle={angle_h:.1f} deg"
            print(line)

    def export_results(self):
        """Export all points and measurements to CSV."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base = os.path.splitext(self.image_path)[0]
        csv_path = f"{base}_measurements_{timestamp}.csv"

        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)

            # Header
            w.writerow(['Citicorp Center — Structural Dimension Extraction'])
            w.writerow([f'Image: {self.image_path}'])
            w.writerow([f'Date: {datetime.now().isoformat()}'])
            if self.ref_scale is not None:
                w.writerow([f'Scale: {self.ref_scale:.4f} px/ft'])
            w.writerow([])

            # Points
            w.writerow(['--- Points ---'])
            w.writerow(['Label', 'X (px)', 'Y (px)', 'X (ft)', 'Y (ft)'])
            for i, (px, py) in enumerate(self.points):
                label = self.labels[i] if i < len(self.labels) else f"P{i+1}"
                if self.ref_scale is not None:
                    w.writerow([label, px, py, f"{px/self.ref_scale:.2f}", f"{py/self.ref_scale:.2f}"])
                else:
                    w.writerow([label, px, py, '', ''])

            w.writerow([])

            # Measurements
            w.writerow(['--- Measurements ---'])
            for m in self.measurements:
                if m['type'] == 'distance':
                    ft_str = f"{m.get('ft', 'N/A'):.2f}" if 'ft' in m else 'N/A'
                    w.writerow(['Distance', m['labels'][0], m['labels'][1],
                               f"{m['px']:.1f} px", f"{ft_str} ft"])
                elif m['type'] == 'angle':
                    w.writerow(['Angle', m['labels'][0], m['labels'][1], m['labels'][2],
                               f"{m['angle_deg']:.1f} deg"])
                elif m['type'] == 'reference':
                    w.writerow(['Reference', f"{m['px']:.1f} px", f"{m['ft']:.1f} ft",
                               f"Scale: {m['scale']:.4f} px/ft"])

        print(f"  Exported to: {csv_path}")

    def run(self):
        """Main event loop."""
        print("\n" + "=" * 60)
        print("  CITICORP CENTER — STRUCTURAL DIMENSION EXTRACTION")
        print("  OpenCV Point Picker")
        print("=" * 60)
        print("\nInstructions:")
        print("  1. First set a REFERENCE: click 2 points with known")
        print("     distance (e.g., building width = 157 ft), then press 'r'")
        print("  2. Click structural points (bracing nodes, stilt bases, etc.)")
        print("  3. Press 'd' to measure distance between last 2 points")
        print("  4. Press 'a' to measure angle at middle of last 3 points")
        print("  5. Press 'e' to export all measurements to CSV")
        print()

        self._redraw()

        while True:
            key = cv2.waitKey(0) & 0xFF

            if key == ord('q') or key == 27:  # q or ESC
                break
            elif key == ord('r'):
                self.set_reference()
            elif key == ord('d'):
                self.measure_distance()
            elif key == ord('a'):
                self.measure_angle()
            elif key == ord('m'):
                self.measure_all_pairs()
            elif key == ord('e'):
                self.export_results()
            elif key == ord('c'):
                self.points.clear()
                self.labels.clear()
                self.measurements.clear()
                print("  All points cleared.")
                self._redraw()
            elif key == ord('z'):
                self.zoom_active = not self.zoom_active
                if self.zoom_active and self.points:
                    self.zoom_center = self.points[-1]
                    print(f"  Zoom ON at {self.zoom_center}")
                else:
                    self.zoom_center = None
                    print("  Zoom OFF")
                self._redraw()

        cv2.destroyAllWindows()

        # Print final summary
        if self.measurements:
            print("\n" + "=" * 60)
            print("  MEASUREMENT SUMMARY")
            print("=" * 60)
            for m in self.measurements:
                if m['type'] == 'distance':
                    ft_str = f" = {m['ft']:.2f} ft" if 'ft' in m else ""
                    print(f"  Distance {m['labels'][0]}->{m['labels'][1]}: "
                          f"{m['px']:.1f} px{ft_str}")
                elif m['type'] == 'angle':
                    print(f"  Angle at {m['labels'][1]}: {m['angle_deg']:.1f} deg")


def main():
    if len(sys.argv) > 1:
        image_path = sys.argv[1]
    else:
        # Default: Citicorp construction photo
        image_path = os.path.join(os.path.dirname(__file__),
                                  "citicorp-construction.png")

    if not os.path.exists(image_path):
        print(f"Error: Image not found: {image_path}")
        print("Usage: python extract_dimensions.py [image_path]")
        sys.exit(1)

    picker = PointPicker(image_path)
    picker.run()


if __name__ == "__main__":
    main()
