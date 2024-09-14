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
function createAxisLine(start, end, color, isLeftGrid) {
  const material = new THREE.LineBasicMaterial({ color: color, depthTest: false });
  const points = [start, end];
  const geometry = new THREE.BufferGeometry().setFromPoints(points);
  const line = new THREE.Line(geometry, material);

  // Инвертируем координаты X для линии оси X левой сетки
  if (isLeftGrid) {
    line.geometry.attributes.position.array[0] *= -1;
    line.geometry.attributes.position.array[3] *= -1;
    line.geometry.attributes.position.needsUpdate = true;
  }

  return line;
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

  let xAxisOffset = isLeftGrid ? gridWidth : 0;

  // --- Группа для оси X ---
  const xAxisGroup = new THREE.Group();
  axisGroup.add(xAxisGroup);

  const xAxis = createSphere(xColor, sphereRadius);
  xAxis.position.set(length, 0, 0);

  if (isLeftGrid) {
    xAxis.position.x *= -1; // Отражаем позицию синей точки для левой сетки
  }

  xAxisGroup.add(xAxis);

  const xAxisLine = createAxisLine(new THREE.Vector3(0, 0, 0), new THREE.Vector3(length, 0, 0), xColor, isLeftGrid);
  xAxisGroup.add(xAxisLine);

  xAxisGroup.position.set(xAxisOffset, 0, 0);

  // --- Группа для оси Y ---
  const yAxisGroup = new THREE.Group();
  axisGroup.add(yAxisGroup);

  const yAxis = createSphere(yColor, sphereRadius);
  yAxis.position.set(0, length, 0); 
  yAxisGroup.add(yAxis);

  const yAxisLine = createAxisLine(new THREE.Vector3(0, 0, 0), new THREE.Vector3(0, length, 0), yColor, isLeftGrid);
  yAxisGroup.add(yAxisLine);

  yAxisGroup.position.set(xAxisOffset, 0, 0);

  // --- Группа для оси Z ---
  const zAxisGroup = new THREE.Group();
  axisGroup.add(zAxisGroup);

  const zAxis = createSphere(zColor, sphereRadius);
  zAxis.position.set(0, 0, length);
  zAxisGroup.add(zAxis);

  const zAxisLine = createAxisLine(new THREE.Vector3(0, 0, 0), new THREE.Vector3(0, 0, length), zColor, isLeftGrid);
  zAxisGroup.add(zAxisLine);

  zAxisGroup.position.set(xAxisOffset, 0, 0);

  // Нулевая точка
  const zeroPoint = createZeroPoint(purpleMaterial);
  zeroPoint.position.set(xAxisOffset, 0, 0); 
  axisGroup.add(zeroPoint); 

  return axisGroup;
}

// --- Создание сеток ---
var grid1 = createGrid(gridWidth, gridHeight, gridDepth, cellSize, 0xffffff);
var grid2 = createGrid(gridWidth, gridHeight, gridDepth, cellSize, 0xffffff);

// --- Создаем общий родительский объект для сеток ---
const sequencerGroup = new THREE.Group();

// --- Создаем родительский объект для левой сетки ---
const leftSequencerGroup = new THREE.Group();

// --- Создаем родительский объект для правой сетки ---
const rightSequencerGroup = new THREE.Group();

// --- Вычисляем координаты центральной точки между фиолетовыми точками ---
const centerX = (leftSequencerGroup.position.x + rightSequencerGroup.position.x) / 2;

// --- Устанавливаем центр вращения ---
sequencerGroup.position.set(centerX, 0, 0);

scene.add(sequencerGroup);

sequencerGroup.add(leftSequencerGroup); // Добавляем в общий родительский объект
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

// ---  Позиционируем левую сетку ---
leftSequencerGroup.position.set(-gridWidth / 2 - gap / 2, -gridHeight / 2, -gridDepth / 2);

// ---  Позиционируем правую сетку ---
rightSequencerGroup.position.set(gridWidth / 2 + gap / 2, -gridHeight / 2, -gridDepth / 2);

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
    mesh.rotation.y = rotation;
    geometry.computeBoundingBox();
    var labelWidth = geometry.boundingBox.max.x - geometry.boundingBox.min.x;
    mesh.position.x -= labelWidth / 2;

    // Добавляем метку в соответствующую группу
    if (isLeftGrid) {
      leftSequencerGroup.add(mesh);
    } else {
      rightSequencerGroup.add(mesh);
    }

    return mesh;
  }

  // --- Создание меток для осей ---
  createLabel('-180°', cellSize * 16, new THREE.Vector3(leftAxis.children[0].position.x, leftAxis.children[0].position.y, leftAxis.children[0].position.z), 0, labelMaterial, true);
  createLabel('+180°', cellSize * 16, new THREE.Vector3(rightAxis.children[0].position.x, rightAxis.children[0].position.y, rightAxis.children[0].position.z), 0, labelMaterial, false);

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

    const intersects = raycaster.intersectObjects([leftAxis.children[0].children[0], leftAxis.children[1].children[0], leftAxis.children[2].children[0], rightAxis.children[0].children[0], rightAxis.children[1].children[0], rightAxis.children[2].children[0]]);

    if (intersects.length > 0) {
      const selectedSphere = intersects[0].object;

      let isLeftGrid = false;

      if (selectedSphere === leftAxis.children[0].children[0]) {
        selectedX = -axisLength;
        selectedY = 0;
        selectedZ = 0;
        isLeftGrid = true;
      } else if (selectedSphere === leftAxis.children[1].children[0]) {
        selectedX = 0;
        selectedY = axisLength;
        selectedZ = 0;
        isLeftGrid = true;
      } else if (selectedSphere === leftAxis.children[2].children[0]) {
        selectedX = 0;
        selectedY = 0;
        selectedZ = axisLength;
        isLeftGrid = true;
      } else if (selectedSphere === rightAxis.children[0].children[0]) {
        selectedX = axisLength;
        selectedY = 0;
        selectedZ = 0;
      } else if (selectedSphere === rightAxis.children[1].children[0]) {
        selectedX = 0;
        selectedY = axisLength;
        selectedZ = 0;
      } else if (selectedSphere === rightAxis.children[2].children[0]) {
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
      currentColumn = createColumn(selectedX, selectedY, selectedZ, greenMaterial, isLeftGrid);

      // --- Добавляем столбец в соответствующую группу ---
      if (isLeftGrid) {
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

    scaleFactor += delta * 0.005; 
    scaleFactor = Math.max(0.25, Math.min(scaleFactor, 2));

    // Масштабируем общий родительский объект
    sequencerGroup.scale.set(scaleFactor, scaleFactor, scaleFactor);

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
        startButton.textContent = "Start"; // Возвращаем текст кнопки на "Start”
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
// --- Функция для создания линии ---
function createLine(x1, y1, z1, x2, y2, z2, color, opacity, isLeftGrid, width, offset) { 
  // Проверка границ индекса
  if (y1 < 1 || y1 > semitones.length || y2 < 1 || y2 > semitones.length) {
    return new THREE.Line(); // Возвращаем пустую линию, если индекс выходит за границы
  }

  if (isLeftGrid) {
    // Отсчитываем координаты X от правой грани (gridWidth) и двигаемся влево
    x1 = gridWidth - x1;
    x2 = gridWidth - x2;

    // Проверка, чтобы левая граница столбца не выходила за пределы сетки
    if (x1 < 0) {
      x1 = 0;
    }
    if (x2 < 0) {
      x2 = 0;
    }

    // Проверка смещения влево и закрашивание синим
    if (offset > 0 && x1 < gridWidth - degreesToCells(semitones[y1 - 1].deg) && x1 >= gridWidth - degreesToCells(semitones[y1 - 1].deg) - offset) {
      color = blueMaterial;
    } else if (x1 >= gridWidth - degreesToCells(semitones[y1 - 1].deg)) { // Внутри зеленой призмы
      color = greenMaterial;
    }
  } else {
    // Проверка смещения вправо и закрашивание красным
    if (offset > 0 && x2 > degreesToCells(semitones[y2 - 1].deg) && x2 <= degreesToCells(semitones[y2 - 1].deg) + offset) {
      color = redMaterial;
    } else if (x2 <= degreesToCells(semitones[y2 - 1].deg)) { // Внутри зеленой призмы
      color = greenMaterial;
    }
  }
  const material = new THREE.LineBasicMaterial({ color: color, opacity: opacity, transparent: true });
  const points = [new THREE.Vector3(x1, y1, z1), new THREE.Vector3(x2, y2, z2)];
  const geometry = new THREE.BufferGeometry().setFromPoints(points);
  const line = new THREE.Line(geometry, material);
  return line;
}

// --- Функция для расчета яркости и прозрачности ---
function calculateOpacityAndColor(dB) {
  const minOpacity = 0.01;
  const maxOpacity = 0.5;
  const minColor = new THREE.Color(0xffffff);
  const maxColor = new THREE.Color(0x00ff00);

  const opacity = minOpacity + (maxOpacity - minOpacity) * (dB / 130);
  const color = minColor.lerp(maxColor, dB / 130);

  return { opacity, color };
}

// --- Функция для создания столбца ---
function createColumn(x, y, z, color, isLeftGrid) {
  const { opacity, color: lineColor } = calculateOpacityAndColor(z); // Используем z как dB

  const width = degreesToCells(semitones[y - 1].deg); // y - 1, так как индексация полутонов начинается с 0
  // Изменено: устанавливаем начальное положение столбца в зависимости от сетки
  const startX = isLeftGrid ? leftAxis.children[2].position.x : 0; // Координата X белой оси
  
  // Создаем группу для линий столбца
  const columnGroup = new THREE.Group();
  columnGroup.position.x = startX; // Устанавливаем начальное положение столбца

  // Создаем линии вокруг ячейки
  const lines = [
    createLine(0, y, 1, width, y, 1, lineColor, opacity, isLeftGrid, width, x),
    createLine(0, y, z, width, y, z, lineColor, opacity, isLeftGrid, width, x),
    createLine(0, y, 1, 0, y, z, lineColor, opacity, isLeftGrid, width, x),
    createLine(width, y, 1, width, y, z, lineColor, opacity, isLeftGrid, width, x),
    createLine(0, y + 1, 1, width, y + 1, 1, lineColor, opacity, isLeftGrid, width, x),
    createLine(0, y + 1, z, width, y + 1, z, lineColor, opacity, isLeftGrid, width, x),
    createLine(0, y + 1, 1, 0, y + 1, z, lineColor, opacity, isLeftGrid, width, x),
    createLine(width, y + 1, 1, width, y + 1, z, lineColor, opacity, isLeftGrid, width, x),
    createLine(0, y, 1, 0, y + 1, 1, lineColor, opacity, isLeftGrid, width, x),
    createLine(width, y, 1, width, y + 1, 1, lineColor, opacity, isLeftGrid, width, x),
    createLine(0, y, z, 0, y + 1, z, lineColor, opacity, isLeftGrid, width, x),
    createLine(width, y, z, width, y + 1, z, lineColor, opacity, isLeftGrid, width, x),
  ];

  // Добавляем линии в группу
  lines.forEach(line => columnGroup.add(line));

  // Возвращаем группу линий (столбец)
  return columnGroup;
}

// --- Преобразование градусов в количество ячеек ---
function degreesToCells(degrees) {
  return Math.round(degrees / 180 * 130);
}

const columns = []; // Массив для хранения столбцов

// --- Отрисовка столбцов для каждого полутона ---
for (let i = 1; i < semitones.length; i++) {
  const randomDB = Math.round(Math.random() * 130);
  // Вычисляем максимальное смещение для текущего столбца
  let maxOffset = gridWidth - degreesToCells(semitones[i].deg);

  let columnLeft = createColumn(0, i + 1, randomDB, greenMaterial, true); // Левая сетка
  let columnRight = createColumn(0, i + 1, randomDB, greenMaterial, false); // Правая сетка
  columns.push({
    left: columnLeft,
    right: columnRight,
    direction: -1,                   // Начальное направление для левой сетки (-1 - влево)
    maxOffset: maxOffset,
    width: degreesToCells(semitones[i].deg),  // Добавлено: ширина столбца
    speed: Math.random() * 0.2 + 0.1, // Скорость движения
    dB: randomDB,                   // Текущий уровень громкости
    dBDirection: 1                   // Направление изменения громкости (1 - вверх, -1 - вниз)
  });
  leftSequencerGroup.add(columnLeft);
  rightSequencerGroup.add(columnRight);
}

// --- Функция для анимации столбцов и рендеринга ---
function animate() {
  requestAnimationFrame(animate);

  // Анимация столбцов
  columns.forEach(columnData => {
    // Движение вдоль оси X
    if (columnData.left.parent) {
      columnData.left.position.x -= columnData.direction * columnData.speed;
      columnData.left.scale.z = columnData.dB / 130;

      // Изменение направления движения для левой сетки
      if (columnData.left.position.x <= 0 || columnData.left.position.x >= columnData.maxOffset) { 
        columnData.direction *= -1; 
      }
    }

    if (columnData.right.parent) {
      columnData.right.position.x += columnData.direction * columnData.speed;
      columnData.right.scale.z = columnData.dB / 130;

      // Изменение направления движения для правой сетки
      if (columnData.right.position.x >= columnData.maxOffset || columnData.right.position.x <= 0) {
        columnData.direction *= -1; 
      }
    }

    // Изменение громкости (dB)
    columnData.dB += columnData.dBDirection * (Math.random() * 0.5 + 0.1);
    if (columnData.dB >= 130 || columnData.dB <= 0) {
      columnData.dBDirection *= -1; 
    }
  });

  // Рендеринг сцены
  renderer.render(scene, camera);
}
animate();