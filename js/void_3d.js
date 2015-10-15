// Void Evolution in Pure Aluminum
// 3D display related functions which based on "Three.js"
// Tianyu Liu

var camera, controls, scene, renderer;

var sGeo = new THREE.SphereGeometry(1.4);
var mpg = new THREE.MeshPhongMaterial({
  color : 0xb2b2b2,
  shading : THREE.FlatShading
});

var x_label = new THREE.Sprite(new THREE.SpriteMaterial({
  map : THREE.ImageUtils.loadTexture("pics/002_64.png")
}));
var y_label = new THREE.Sprite(new THREE.SpriteMaterial({
  map : THREE.ImageUtils.loadTexture("pics/020_64.png")
}));
var z_label = new THREE.Sprite(new THREE.SpriteMaterial({
  map : THREE.ImageUtils.loadTexture("pics/200_64.png")
}));
var lMat = new THREE.LineBasicMaterial()
var cMesh = new THREE.Mesh(new THREE.CubeGeometry(a, a, a), lMat);
unitCell = new THREE.BoxHelper(cMesh);
unitCell.material.color.set(0xffffff);

var bAtom = new THREE.PointsMaterial({
  size : 30,
  map : THREE.ImageUtils.loadTexture("pics/blue_atom_64.png"),
  alphaTest : 0.5,
  sizeAttenuation : false
});

var light = new THREE.AmbientLight(0xffffff);
// var light = new THREE.DirectionalLight(0xffffff);
// light.position.set(0, 0, 1);

function init() {

  camera = new THREE.OrthographicCamera(-a, a, a, -a, -a * 2, 2 * a);
  camera.position.z = a;

  controls = new THREE.OrthographicTrackballControls(camera, $("container"));
  controls.zoomSpeed = 10;
  controls.staticMoving = false;
  controls.dynamicDampingFactor = 1;
  controls.addEventListener('change', render);
  // controls.addEventListener('change', light_update);

  // world
  scene = new THREE.Scene();
  scene.add(light);
  scene.add(unitCell);

  var axisHelper = new THREE.AxisHelper(a);
  scene.add(axisHelper);

  // renderer
  renderer = new THREE.WebGLRenderer();
  renderer.setPixelRatio(window.devicePixelRatio);
  renderer.setSize(400, 400);
  $("3d").appendChild(renderer.domElement);

  // compass
  renderer2 = new THREE.WebGLRenderer({
    alpha : true
  });
  renderer2.setSize(100, 100);
  $("compass").appendChild(renderer2.domElement);
  scene2 = new THREE.Scene();
  camera2 = new THREE.OrthographicCamera(-2, 2, 2, -2, -4, 4);
  camera2.up = camera.up;
  scene2.add(camera2);
  var axisHelper = new THREE.AxisHelper(1.5);
  scene2.add(axisHelper);
  x_label.position.set(0, 0, 1.5);
  y_label.position.set(0, 1.5, 0);
  z_label.position.set(1.5, 0, 0);
  scene2.add(x_label);
  scene2.add(y_label);
  scene2.add(z_label);

  animate();
  render();
}

function draw3D(P) {
  var D3 = parseInt($("3d_opt").value), Cd = P[0], CdN = P[1], mL = P[2], mL1 =
      mL * a * 2, i = 0;

  mL0 = mL1 / mL0;

  controls.left0 = camera.left *= mL0;
  controls.bottom0 = camera.bottom *= mL0;
  controls.right0 = camera.right *= mL0;
  controls.top0 = camera.top *= mL0;
  camera.near *= mL0, camera.far *= mL0;

  camera.position.x *= mL0, camera.position.y *= mL0, camera.position.z *= mL0;
  controls.position0.set(0, 0, mL1);
  camera.updateProjectionMatrix();

  n_D3 = D3.length - 1;

  var pGeo = new THREE.Geometry();

  if (D3 === 0) {
    pGeo.vertices = Array(CdN[0]);
    for (i = 0; i < CdN[0]; i++) {
      var pVec3 = new THREE.Vector3;
      pVec3.x = Cd[0][i][0] * a, pVec3.y = Cd[0][i][1] * a, pVec3.z =
          Cd[0][i][2] * a;
      pGeo.vertices[i] = pVec3;
    }
  } else {
    var n = 0, ii = 0;
    for (i = D3; i < 16; i++)
      n += CdN[i];

    pGeo.vertices = Array(n);
    for (j = D3; j < 16; j++)
      for (i = 0; i < CdN[j]; i++) {
        // console.log(Cd[j][i],j,i, CdN[j])
        var pVec3 = new THREE.Vector3;
        pVec3.x = Cd[j][i][0] * a, pVec3.y = Cd[j][i][1] * a, pVec3.z =
            Cd[j][i][2] * a;
        pGeo.vertices[ii++] = pVec3;
      }
  }

  scene.children[1].scale.set(mL * 2, mL * 2, mL * 2);
  scene.children[2] = new THREE.AxisHelper(mL1);
  bAtom.size = 600 / mL1 * camera.zoom;
  scene.children[3] = new THREE.Points(pGeo, bAtom);

  mL0 = mL1;

  animate()
  render();
}

function light_update() {
  light.position.copy(camera.position);
}

function zoom(scale) {
  (camera.zoom) *= scale;
  camera.updateProjectionMatrix();
  if (scene.children[3] != null) scene.children[3].material.size *= scale;

  render();
}

function resetCam() {
  controls.reset();
  animate()
  render()
}

function animate() {
  requestAnimationFrame(animate);
  controls.update();
  camera2.position.x = (camera.position.x - controls.target.x) / mL0 * 2;
  camera2.position.y = (camera.position.y - controls.target.y) / mL0 * 2;
  camera2.position.z = (camera.position.z - controls.target.z) / mL0 * 2;
  camera2.position.setLength(2);
  camera2.lookAt(scene2.position);
}

function render() {
  renderer.render(scene, camera);
  renderer2.render(scene2, camera2);
}