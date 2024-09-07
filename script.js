// --- Размеры сетки (130x130x130) ---
const gridWidth = 130;
const gridHeight = 130;
const gridDepth = 130;

// --- Размер ячейки ---
const cellSize = 1;

// --- Промежуток между сетками ---
const gap = 5;

// --- Three.js setup ---
const scene = new THREE.Scene();
const camera = new THREE.OrthographicCamera(-window.innerWidth / 2, window.innerWidth / 2, window.innerHeight / 2, -window.innerHeight / 2, -10000, 10000);
const renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: "high-performance" });
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(window.innerWidth, window.innerHeight);
document.getElementById('grid-container').appendChild(renderer.domElement);

// --- Материалы ---
const defaultMaterial = new THREE.LineBasicMaterial({
  color: 0xffffff,
  opacity: 0.01, // Уменьшаем непрозрачность
  transparent: true,
  depthWrite: false,
  depthTest: false
});

const greenMaterial = new THREE.Color(0x00ff00);
const blueMaterial = new THREE.Color(0x0000ff);
const redMaterial = new THREE.Color(0xff0000);

// --- Кислотно-фиолетовый цвет ---
const purpleMaterial = new THREE.Color(0x800080);

const labelMaterial = new THREE.MeshBasicMaterial({ color: 0xffffff, depthTest: false });

// --- Функция для создания шарика ---
function createSphere(color, radius) {
  const geometry = new THREE.SphereGeometry(radius, 32, 32);
  const material = new THREE.MeshBasicMaterial({ color: color, depthTest: false });
  const sphere = new THREE.Mesh(geometry, material);
  return sphere;
}

// --- Функция для создания нулевой точки ---
function createZeroPoint(color) {
  const geometry = new THREE.SphereGeometry(sphereRadius, 32, 32);
  const material = new THREE.MeshBasicMaterial({ color: color, depthTest: false });
  const sphere = new THREE.Mesh(geometry, material);
  return sphere;
}

// --- Функция для создания линии оси ---
function createAxisLine(start, end, color) {
  const material = new THREE.LineBasicMaterial({ color: color, depthTest: false });
  const points = [start, end];
  const geometry = new THREE.BufferGeometry().setFromPoints(points);
  const line = new THREE.Line(geometry, material);
  return line;
}

// --- Функция для создания столбца ---
function createColumn(x, y, z, color, isLeftGrid) {
  const geometry = new THREE.BoxGeometry(cellSize, cellSize, z * cellSize);
  const material = new THREE.MeshBasicMaterial({ color: color });
  const column = new THREE.Mesh(geometry, material);
  // Позиционируем столбец относительно родительской группы
  if (isLeftGrid) {
    column.position.set(x + leftSequencerGroup.position.x, y + leftSequencerGroup.position.y, z / 2 + leftSequencerGroup.position.z);
  } else {
    column.position.set(x + rightSequencerGroup.position.x, y + rightSequencerGroup.position.y, z / 2 + rightSequencerGroup.position.z);
  }
  return column;
}

// --- Функция для создания сетки ---
function createGrid(gridWidth, gridHeight, gridDepth, cellSize, color) {
  const geometry = new THREE.BufferGeometry();
  const positions = [];

  // --- Создание линий сетки ---
  for (let y = 0; y <= gridHeight; y += 1) {
    for (let z = 0; z <= gridDepth; z += 1) {
      positions.push(0, y * cellSize, z * cellSize, gridWidth * cellSize, y * cellSize, z * cellSize);
    }
  }

  for (let x = 0; x <= gridWidth; x += 1) {
    for (let z = 0; z <= gridDepth; z += 1) {
      positions.push(x * cellSize, 0, z * cellSize, x * cellSize, gridHeight * cellSize, z * cellSize);
    }
  }

  for (let z = 0; z <= gridWidth; z += 1) {
    for (let y = 0; y <= gridHeight; y += 1) {
      positions.push(z * cellSize, y * cellSize, 0, z * cellSize, y * cellSize, gridDepth * cellSize);
    }
  }

  geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));

  var material = new THREE.LineBasicMaterial({ color: color, opacity: 0.01, transparent: true, depthWrite: false, depthTest: false }); // Уменьшаем непрозрачность
  return new THREE.LineSegments(geometry, material);
}

// --- Функция для создания осей ---
function createAxis(length, sphereRadius, xColor, yColor, zColor, isLeftGrid) {
  const axisGroup = new THREE.Group();

  const xAxis = createSphere(xColor, sphereRadius);
  xAxis.position.set(length, 0, 0);
  axisGroup.add(xAxis);

  const yAxis = createSphere(yColor, sphereRadius);
  yAxis.position.set(0, length, 0); 
  axisGroup.add(yAxis);

  const zAxis = createSphere(zColor, sphereRadius);
  zAxis.position.set(0, 0, length);
  axisGroup.add(zAxis);

  const zeroPoint = createZeroPoint(purpleMaterial);
  zeroPoint.position.set(0, 0, 0);
  axisGroup.add(zeroPoint);

  const xAxisLine = createAxisLine(zeroPoint.position, xAxis.position, xColor);
  axisGroup.add(xAxisLine);

  const yAxisLine = createAxisLine(zeroPoint.position, yAxis.position, yColor);
  axisGroup.add(yAxisLine);

  const zAxisLine = createAxisLine(zeroPoint.position, zAxis.position, zColor);
  axisGroup.add(zAxisLine);

  // Отражаем позицию точки X для левой сетки
  if (isLeftGrid) {
    xAxis.position.x *= -1;
  }

  return axisGroup;
}
// --- Создание сеток ---
var grid1 = createGrid(gridWidth, gridHeight, gridDepth, cellSize, 0xffffff);
var grid2 = createGrid(gridWidth, gridHeight, gridDepth, cellSize, 0xffffff);

// --- Создаем общий родительский объект для сеток ---
const sequencerGroup = new THREE.Group();
scene.add(sequencerGroup);

// --- Создаем родительский объект для левой сетки ---
const leftSequencerGroup = new THREE.Group();
sequencerGroup.add(leftSequencerGroup); // Добавляем в общий родительский объект

// --- Создаем родительский объект для правой сетки ---
const rightSequencerGroup = new THREE.Group();
sequencerGroup.add(rightSequencerGroup); // Добавляем в общий родительский объект

// --- Создание осей с шариками ---
var axisLength = 130; // Длина каждой оси 130 ячеек
var sphereRadius = 5; // Радиус шариков

// --- Создаем оси ---
const leftAxis = createAxis(axisLength, sphereRadius, blueMaterial, greenMaterial, 0xffffff, true);
leftSequencerGroup.add(leftAxis);

const rightAxis = createAxis(axisLength, sphereRadius, redMaterial, greenMaterial, 0xffffff, false);
rightSequencerGroup.add(rightAxis);

// --- Добавляем сетки в соответствующие группы ---
leftSequencerGroup.add(grid1);
rightSequencerGroup.add(grid2);

// --- Отражаем левую сетку ---
leftSequencerGroup.scale.set(-1, 1, 1);
leftSequencerGroup.position.set(-gridWidth / 2 - gap / 2, -gridHeight / 2, -gridDepth / 2);

// --- Центрируем правую сетку ---
rightSequencerGroup.position.set(gridWidth / 2 + gap / 2, -gridHeight / 2, -gridDepth / 2);

// --- Центрируем секвенсор для вращения ---
sequencerGroup.position.set(0, 0, 0); // Центрируем группу секвенсора

// --- Переменные для хранения координат выбранной ячейки ---
let selectedX = 0;
let selectedY = 0;
let selectedZ = 0;

// --- Создаем столбец ---
let currentColumn = null;

// --- Добавление меток ---
var loader = new THREE.FontLoader();
loader.load('https://threejs.org/examples/fonts/helvetiker_regular.typeface.json', function (font) {
  // --- Функция для создания меток ---
  function createLabel(text, size, position, rotation = 0, labelMaterial, isLeftGrid) {
    var geometry = new THREE.TextGeometry(text, {
      font: font,
      size: size,
      height: 0.1
    });
    var mesh = new THREE.Mesh(geometry, labelMaterial);
    mesh.position.copy(position);

    // Зеркально отражаем позицию текста для левой сетки
    if (isLeftGrid) {
      mesh.position.x *= -1;
      geometry.computeBoundingBox();
      var labelWidth = geometry.boundingBox.max.x - geometry.boundingBox.min.x;
      mesh.position.x -= labelWidth; // Корректируем позицию после отражения
    }
    mesh.rotation.y = rotation;

    // Добавляем метку в соответствующую группу
    if (isLeftGrid) {
      leftSequencerGroup.add(mesh);
    } else {
      rightSequencerGroup.add(mesh);
    }

    return mesh;
  }

  // --- Создание меток для осей ---
  createLabel('-180°', cellSize * 16, new THREE.Vector3(leftAxis.children[0].position.x - cellSize * 15, leftAxis.children[0].position.y, leftAxis.children[0].position.z), 0, labelMaterial, true);
  createLabel('+180°', cellSize * 16, new THREE.Vector3(rightAxis.children[0].position.x + cellSize * 15, rightAxis.children[0].position.y, rightAxis.children[0].position.z), 0, labelMaterial, false);

  createLabel(semitones[0].f.toFixed(2), cellSize * 8, new THREE.Vector3(leftAxis.children[2].position.x, leftAxis.children[2].position.y, leftAxis.children[2].position.z + cellSize * 10), 0, labelMaterial, true);
  createLabel(semitones[semitones.length - 1].f.toFixed(2), cellSize * 8, new THREE.Vector3(rightAxis.children[2].position.x, rightAxis.children[2].position.y, rightAxis.children[2].position.z + cellSize * 10), 0, labelMaterial, false);

  createLabel('0', cellSize * 8, new THREE.Vector3(leftAxis.children[3].position.x + cellSize * 5, leftAxis.children[3].position.y, leftAxis.children[3].position.z), 0, labelMaterial, true);
  createLabel('0', cellSize * 8, new THREE.Vector3(rightAxis.children[3].position.x + cellSize * 5, rightAxis.children[3].position.y, rightAxis.children[3].position.z), 0, labelMaterial, false);

  // --- Обработчик события touchstart для выбора осей ---
  renderer.domElement.addEventListener('touchstart', function (event) {
    event.preventDefault();

    const rect = renderer.domElement.getBoundingClientRect();
    const mouse = new THREE.Vector2();
    mouse.x = ((event.touches[0].pageX - rect.left) / renderer.domElement.clientWidth) * 2 - 1;
    mouse.y = -((event.touches[0].pageY - rect.top) / renderer.domElement.clientHeight) * 2 + 1;

    const raycaster = new THREE.Raycaster();
    raycaster.setFromCamera(mouse, camera);

    const intersects = raycaster.intersectObjects([leftAxis.children[0], leftAxis.children[1], leftAxis.children[2], rightAxis.children[0], rightAxis.children[1], rightAxis.children[2]]);

    if (intersects.length > 0) {
      const selectedSphere = intersects[0].object;

      if (selectedSphere === leftAxis.children[0]) {
        selectedX = -axisLength;
        selectedY = 0;
        selectedZ = 0;
      } else if (selectedSphere === leftAxis.children[1]) {
        selectedX = 0;
        selectedY = axisLength;
        selectedZ = 0;
      } else if (selectedSphere === leftAxis.children[2]) {
        selectedX = 0;
        selectedY = 0;
        selectedZ = axisLength;
      } else if (selectedSphere === rightAxis.children[0]) {
        selectedX = axisLength;
        selectedY = 0;
        selectedZ = 0;
      } else if (selectedSphere === rightAxis.children[1]) {
        selectedX = 0;
        selectedY = axisLength;
        selectedZ = 0;
      } else if (selectedSphere === rightAxis.children[2]) {
        selectedX = 0;
        selectedY = 0;
        selectedZ = axisLength;
      }

      // --- Удаляем предыдущий столбец ---
      if (currentColumn) {
        if (currentColumn.parent === leftSequencerGroup) {
          leftSequencerGroup.remove(currentColumn);
        } else {
          rightSequencerGroup.remove(currentColumn);
        }
      }
      // --- Создаем новый столбец ---
      currentColumn = createColumn(selectedX, selectedY, selectedZ, greenMaterial, selectedSphere === leftAxis.children[0] || selectedSphere === leftAxis.children[1] || selectedSphere === leftAxis.children[2]);

      // --- Добавляем столбец в соответствующую группу ---
      if (selectedSphere === leftAxis.children[0] || selectedSphere === leftAxis.children[1] || selectedSphere === leftAxis.children[2]) {
        leftSequencerGroup.add(currentColumn);
      } else {
        rightSequencerGroup.add(currentColumn);
      }
    }
  });
});

// --- Установка начального положения камеры ---
camera.position.set(0, 0, 2000);
camera.lookAt(0, 0, 0);

// --- Уменьшаем начальный масштаб ---
scaleFactor = 0.6;
leftSequencerGroup.scale.set(scaleFactor, scaleFactor, scaleFactor);
rightSequencerGroup.scale.set(scaleFactor, scaleFactor, scaleFactor);

// --- Вращение и масштабирование с помощью сенсорного управления ---
var previousTouchX = null;
var previousTouchY = null;
var previousDistance = null;
var isDragging = false;

// Scale factor
var scaleFactor = 1;

renderer.domElement.addEventListener('touchstart', function (event) {
  if (event.touches.length === 1) {
    const rect = renderer.domElement.getBoundingClientRect();
    previousTouchX = event.touches[0].pageX - rect.left;
    previousTouchY = event.touches[0].pageY - rect.top;
    isDragging = true;
  } else if (event.touches.length === 2) {
    var dx = event.touches[0].clientX - event.touches[1].clientX;
    var dy = event.touches[0].clientY - event.touches[1].clientY;
    previousDistance = Math.sqrt(dx * dx + dy * dy);
    isDragging = false;
  }
});

renderer.domElement.addEventListener('touchmove', function (event) {
  if (isDragging && event.touches.length === 1 && previousTouchX !== null && previousTouchY !== null) {
    const rect = renderer.domElement.getBoundingClientRect();
    var deltaX = (event.touches[0].pageX - rect.left) - previousTouchX;
    var deltaY = (event.touches[0].pageY - rect.top) - previousTouchY;

    // Ограничение вращения по осям X и Y (45 градусов)
    const rotationLimit = THREE.MathUtils.degToRad(45);

    // Вращаем только по осям X и Y
    sequencerGroup.rotation.y = Math.max(-rotationLimit, Math.min(sequencerGroup.rotation.y + deltaX * 0.01, rotationLimit));
    sequencerGroup.rotation.x = Math.max(-rotationLimit, Math.min(sequencerGroup.rotation.x + deltaY * 0.01, rotationLimit));

    previousTouchX = event.touches[0].pageX - rect.left;
    previousTouchY = event.touches[0].pageY - rect.top;

  } else if (!isDragging && event.touches.length === 2 && previousDistance) {
    var dx = event.touches[0].clientX - event.touches[1].clientX;
    var dy = event.touches[0].clientY - event.touches[1].clientY;
    var distance = Math.sqrt(dx * dx + dy * dy);
    var delta = distance - previousDistance;

    scaleFactor += delta * 0.01;
    scaleFactor = Math.max(0.25, Math.min(scaleFactor, 1));

    // Масштабируем общий родительский объект
    if (event.touches[0].pageX > window.innerWidth / 2 && event.touches[1].pageX > window.innerWidth / 2) {
      sequencerGroup.scale.set(scaleFactor, scaleFactor, scaleFactor);
    } else if (event.touches[0].pageX < window.innerWidth / 2 && event.touches[1].pageX < window.innerWidth / 2) {
      sequencerGroup.scale.set(scaleFactor, scaleFactor, scaleFactor);
    }

    previousDistance = distance;
  }
});

renderer.domElement.addEventListener('touchend', function (event) {
  previousTouchX = null;
  previousTouchY = null;
  previousDistance = null;
  isDragging = false;
});

// --- Удаляем разделительную линию ---

// Цикл рендеринга
function animate() {
  requestAnimationFrame(animate);
  renderer.render(scene, camera);
}
animate();

// --- Обработчик изменения размера окна ---
window.addEventListener('resize', function () {
  camera.left = -window.innerWidth / 2;
  camera.right = window.innerWidth / 2;
  camera.top = window.innerHeight / 2;
  camera.bottom = -window.innerHeight / 2;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);

  // Центрируем секвенсор
  sequencerGroup.position.set(0, -gridHeight / 2, -gridDepth / 2);
});


// --- Функция DWT (дискретное вейвлет-преобразование) ---
function DWT(data) {
  const N = data.length;
  let temp = new Float32Array(N);

  // Количество уровней разложения
  const levels = 8;

  for (let level = 0; level < levels; level++) {
    const half = N >> (level + 1);

    // Вычисление коэффициентов аппроксимации
    for (let i = 0; i < half; i++) {
      temp[i] = (data[2 * i] + data[2 * i + 1]) / 2;
    }

    // Вычисление коэффициентов детализации
    for (let i = 0; i < half; i++) {
      temp[half + i] = (data[2 * i] - data[2 * i + 1]) / 2;
    }

    // Копирование временных данных обратно в исходный массив
    for (let i = 0; i < N; i++) {
      data[i] = temp[i];
    }
  }

  return data;
}

// --- Обработчик события click кнопки "Start/Stop" ---
const startButton = document.getElementById('startButton');
let isRecording = false; // Флаг для отслеживания состояния записи
let audioContext, source, analyser, processAudioInterval; // Объявляем переменные для аудио

startButton.addEventListener('click', function () {
  if (!isRecording) {
    console.log("Start button clicked!");
    startButton.textContent = "Stop"; // Изменяем текст кнопки на "Stop"
    navigator.mediaDevices.getUserMedia({ audio: true })
      .then(stream => {
        audioContext = new AudioContext();
        source = audioContext.createMediaStreamSource(stream);
        analyser = audioContext.createAnalyser();
        source.connect(analyser);

        const dataArray = new Uint8Array(analyser.frequencyBinCount);

        // Опорная амплитуда для расчета dB
        const referenceAmplitude = 1;

        function processAudio() {
          analyser.getByteTimeDomainData(dataArray);

          // Применяем DWT
          const dwtCoefficients = DWT(dataArray);

          const semitoneAmplitudes = [];
          for (let i = 0; i < semitones.length; i++) {
            const semitone = semitones[i];

            // Определение диапазона частот для полутона
            const freqStart = semitone.f;
            const freqEnd = i < semitones.length - 1 ? semitones[i + 1].f : audioContext.sampleRate / 2;

            // Находим соответствующие индексы в dwtCoefficients
            const indexStart = Math.floor(freqStart / (audioContext.sampleRate / 2) * dwtCoefficients.length);
            const indexEnd = Math.floor(freqEnd / (audioContext.sampleRate / 2) * dwtCoefficients.length);

            // Проверка на корректность индексов
            if (indexStart >= 0 && indexEnd <= dwtCoefficients.length && indexEnd > indexStart) {
              let sum = 0;
              for (let j = indexStart; j < indexEnd; j++) {
                sum += Math.abs(dwtCoefficients[j]);
              }
              let averageAmplitude = sum / (indexEnd - indexStart);

              averageAmplitude /= 255;

              const dB = 20 * Math.log10(averageAmplitude / referenceAmplitude);

              // Округляем dB до двух знаков после запятой
              semitoneAmplitudes.push(dB.toFixed(2));
            } else {
              // Если индексы некорректны, добавляем 0 dB
              semitoneAmplitudes.push("0.00");
            }
          }

          console.log(semitoneAmplitudes);
        }

        // Запускаем интервал для processAudio
        processAudioInterval = setInterval(processAudio, 50);
        isRecording = true;
      })
      .catch(err => {
        console.error("Error accessing microphone:", err);
        startButton.textContent = "Start"; // Возвращаем текст кнопки на "Start"
        isRecording = false;
      });
  } else {
    console.log("Stop button clicked!");
    startButton.textContent = "Start"; // Изменяем текст кнопки на "Start"

    // Останавливаем захват звука
    if (audioContext) {
      audioContext.close();
      clearInterval(processAudioInterval);
      isRecording = false;
    }
  }
});