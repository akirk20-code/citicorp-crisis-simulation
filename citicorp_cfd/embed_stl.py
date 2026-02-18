#!/usr/bin/env python3
"""Embed STL data as base64 in the dashboard HTML for self-contained 3D viewing."""
import struct, re, base64, os

def ascii_stl_to_binary(ascii_data):
    """Convert ASCII STL to binary STL bytes."""
    text = ascii_data.decode('ascii', errors='ignore')
    triangles = []
    facet_pattern = re.compile(
        r'facet\s+normal\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\n'
        r'\s*outer\s+loop\s*\n'
        r'\s*vertex\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\n'
        r'\s*vertex\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\n'
        r'\s*vertex\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\n'
        r'\s*endloop\s*\n'
        r'\s*endfacet', re.MULTILINE)
    for m in facet_pattern.finditer(text):
        vals = [float(m.group(i)) for i in range(1, 13)]
        triangles.append(vals)
    buf = bytearray(80)
    buf += struct.pack('<I', len(triangles))
    for tri in triangles:
        buf += struct.pack('<12fH', *tri, 0)
    return bytes(buf)

# Paths
script_dir = os.path.dirname(os.path.abspath(__file__))
stl_dir = os.path.join(script_dir, 'constant', 'triSurface')
html_path = os.path.join(script_dir, 'cfd_dashboard.html')

# Convert STLs to binary and base64
b64 = {}
for fname in ['citicorp_tower.stl', 'citicorp_stilts.stl', 'surroundings.stl']:
    with open(os.path.join(stl_dir, fname), 'rb') as f:
        binary = ascii_stl_to_binary(f.read())
    key = fname.replace('.stl', '').replace('citicorp_', '')
    b64[key] = base64.b64encode(binary).decode()
    print(f'{fname}: {len(binary)} bytes binary -> {len(b64[key])} chars base64')

# Build new viewer script
viewer_js = '''// === 3D STL Viewer (Self-contained with embedded STL data) ===
let scene, camera, renderer, controls;
let meshes = { tower: null, stilts: null, surroundings: null };
let wireframeMode = false;

// Base64-encoded binary STL data
const STL_B64 = {
  tower: "''' + b64['tower'] + '''",
  stilts: "''' + b64['stilts'] + '''",
  surroundings: "''' + b64['surroundings'] + '''"
};

function b64ToArrayBuffer(b64) {
  const bin = atob(b64);
  const buf = new ArrayBuffer(bin.length);
  const view = new Uint8Array(buf);
  for (let i = 0; i < bin.length; i++) view[i] = bin.charCodeAt(i);
  return buf;
}

function initViewer() {
  const container = document.getElementById('stl-viewer');
  if (!container) { console.error('stl-viewer container not found'); return; }

  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x070b14);

  camera = new THREE.PerspectiveCamera(45, container.clientWidth / container.clientHeight, 1, 5000);
  camera.position.set(-400, 300, 400);

  renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(container.clientWidth, container.clientHeight);
  renderer.setPixelRatio(window.devicePixelRatio);
  container.appendChild(renderer.domElement);

  controls = new THREE.OrbitControls(camera, renderer.domElement);
  controls.target.set(0, 100, 0);
  controls.enableDamping = true;
  controls.dampingFactor = 0.1;

  // Lights
  scene.add(new THREE.AmbientLight(0x404060, 0.6));
  const dir1 = new THREE.DirectionalLight(0xffffff, 0.8);
  dir1.position.set(200, 400, 300);
  scene.add(dir1);
  const dir2 = new THREE.DirectionalLight(0x4488cc, 0.3);
  dir2.position.set(-200, 200, -100);
  scene.add(dir2);

  // Ground plane
  const ground = new THREE.Mesh(
    new THREE.PlaneGeometry(2000, 2000),
    new THREE.MeshStandardMaterial({ color: 0x0a1020, roughness: 0.9 })
  );
  ground.rotation.x = -Math.PI / 2;
  ground.position.y = -0.5;
  scene.add(ground);

  // Grid
  const grid = new THREE.GridHelper(1500, 30, 0x1a2540, 0x101828);
  scene.add(grid);

  // Wind arrow
  const arrowDir = new THREE.Vector3(1, 0, 0);
  const arrowOrigin = new THREE.Vector3(-350, 100, 0);
  const arrow = new THREE.ArrowHelper(arrowDir, arrowOrigin, 100, 0xff4444, 15, 10);
  scene.add(arrow);

  // Domain box (wireframe)
  const domainGeo = new THREE.BoxGeometry(720, 450, 720);
  const domainMat = new THREE.MeshBasicMaterial({ color: 0x2a3050, wireframe: true, transparent: true, opacity: 0.15 });
  const domainBox = new THREE.Mesh(domainGeo, domainMat);
  domainBox.position.set(0, 225, 0);
  scene.add(domainBox);

  // Parse embedded STL data
  const loader = new THREE.STLLoader();

  try {
    const towerGeo = loader.parse(b64ToArrayBuffer(STL_B64.tower));
    towerGeo.computeVertexNormals();
    towerGeo.rotateX(-Math.PI / 2);
    meshes.tower = new THREE.Mesh(towerGeo, new THREE.MeshPhongMaterial({ color: 0x3b82f6, shininess: 60, flatShading: true }));
    scene.add(meshes.tower);
    console.log('Tower loaded:', towerGeo.attributes.position.count / 3, 'triangles');
  } catch(e) { console.error('Tower STL parse error:', e); }

  try {
    const stiltsGeo = loader.parse(b64ToArrayBuffer(STL_B64.stilts));
    stiltsGeo.computeVertexNormals();
    stiltsGeo.rotateX(-Math.PI / 2);
    meshes.stilts = new THREE.Mesh(stiltsGeo, new THREE.MeshPhongMaterial({ color: 0x06b6d4, shininess: 60, flatShading: true }));
    scene.add(meshes.stilts);
    console.log('Stilts loaded:', stiltsGeo.attributes.position.count / 3, 'triangles');
  } catch(e) { console.error('Stilts STL parse error:', e); }

  try {
    const surrGeo = loader.parse(b64ToArrayBuffer(STL_B64.surroundings));
    surrGeo.computeVertexNormals();
    surrGeo.rotateX(-Math.PI / 2);
    meshes.surroundings = new THREE.Mesh(surrGeo, new THREE.MeshPhongMaterial({ color: 0x64748b, shininess: 20, flatShading: true, transparent: true, opacity: 0.85 }));
    scene.add(meshes.surroundings);
    console.log('Surroundings loaded:', surrGeo.attributes.position.count / 3, 'triangles');
  } catch(e) { console.error('Surroundings STL parse error:', e); }

  animate();
  window.addEventListener('resize', onResize);
  console.log('3D viewer initialized successfully');
}

function animate() {
  requestAnimationFrame(animate);
  controls.update();
  renderer.render(scene, camera);
}

function onResize() {
  const c = document.getElementById('stl-viewer');
  camera.aspect = c.clientWidth / c.clientHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(c.clientWidth, c.clientHeight);
}

function resetCamera() {
  camera.position.set(-400, 300, 400);
  controls.target.set(0, 100, 0);
  controls.update();
  document.querySelectorAll('.viewer-controls button').forEach(b => b.classList.remove('active'));
  event.target.classList.add('active');
}
function viewTop() {
  camera.position.set(0, 600, 0);
  controls.target.set(0, 0, 0);
  controls.update();
  document.querySelectorAll('.viewer-controls button').forEach(b => b.classList.remove('active'));
  event.target.classList.add('active');
}
function viewFront() {
  camera.position.set(-500, 150, 0);
  controls.target.set(0, 150, 0);
  controls.update();
  document.querySelectorAll('.viewer-controls button').forEach(b => b.classList.remove('active'));
  event.target.classList.add('active');
}
function viewSide() {
  camera.position.set(0, 150, 500);
  controls.target.set(0, 150, 0);
  controls.update();
  document.querySelectorAll('.viewer-controls button').forEach(b => b.classList.remove('active'));
  event.target.classList.add('active');
}
function toggleWireframe() {
  wireframeMode = !wireframeMode;
  Object.values(meshes).forEach(m => { if (m) m.material.wireframe = wireframeMode; });
  const btn = document.querySelector('.viewer-controls button:last-child');
  if (wireframeMode) btn.classList.add('active'); else btn.classList.remove('active');
}

initViewer();'''

# Read HTML, replace viewer section
with open(html_path, 'r', encoding='utf-8') as f:
    html = f.read()

old_start = html.find('// === 3D STL Viewer ===')
old_end = html.find('initViewer();', old_start) + len('initViewer();')

if old_start == -1:
    print('ERROR: Could not find old viewer script section')
else:
    new_html = html[:old_start] + viewer_js + html[old_end:]
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(new_html)
    print(f'\nDashboard updated: {len(new_html):,} bytes')
    print('3D viewer now self-contained - no file loading needed!')
    print('Works from file:// or any HTTP server.')
    print('Buttons (Perspective/Top/Front/Side/Wireframe) now update active state.')
