# Профилирование для таблицы 5

## Сохранённый эксперимент и добавление GPU

В ветке `occupancy_test` сохранены исходники, патч HeavyRT, десять JSON-записей
и журналы запуска RX 6800 XT в `profiling/experiments/rx6800xt/`. Для Intel Arc
A380 записи и журналы лежат в `profiling/experiments/arc-a380/`. Общий отчёт
сохранён в `profiling/reports/supplement.tex` и `supplement-document.pdf`.
В начале отчёта описано используемое Vulkan-расширение. Сводная таблица
максимумов по ядрам показывает только счётчики векторных регистров, со столбцами `GPU/PT`
и `GPU/AO`. Вторая сводка со столбцами `GPU/PT` показывает варианты
трассировки со своим именем ядра. Для wavefront это отдельное ядро
`kernel_RayTrace2`, для megakernel — весь объединённый шейдер `PathTraceMega`
вместе с трассировкой и шейдингом. RADV показывает число выделенных VGPR,
а ANV — пиковое давление на GRF до аллокации и spills. Смысл счётчиков
объясняется в подписях к таблицам. Остальные метрики есть в детальных таблицах.

После подготовки зависимостей, описанной ниже, на новой машине:

```bash
git switch occupancy_test
python3 testing/profile_gpu.py build --slicer /path/to/kernel_slicer --rebuild
python3 testing/profile_gpu.py run --integrator mispt --gpu-label RTX4090 \
  --scene test035=scenes/test_035/statex_00001.xml \
  --output profiling_results/rtx4090/test035-pt
python3 testing/profile_gpu.py run --integrator rtao --gpu-label RTX4090 \
  --scene test035=scenes/test_035/statex_00001.xml \
  --output profiling_results/rtx4090/test035-rtao
```

Чтобы добавить результаты в Git без бинарников, изображений и промежуточных
артефактов, скопируйте из каждого каталога запусков файлы `*.json` **кроме**
`*.raw.json` в соответствующие подкаталоги
`profiling/experiments/rtx4090/test035-pt/` и `test035-rtao/`.
Обычная JSON-запись уже содержит полные raw statistics. Журналы `*.log` также
можно сохранить. Затем пересоздайте общий отчёт:

```bash
python3 testing/profile_gpu.py report profiling/experiments \
  --output profiling/reports/supplement.tex
cd profiling/reports
pdflatex -halt-on-error supplement-document.tex
pdflatex -halt-on-error supplement-document.tex
```

`report` читает каталоги рекурсивно: новые GPU становятся отдельными столбцами.
GPU/driver/host различаются автоматически; `--gpu-label` добавляет удобную метку.
Эксперименты с различными сценами, хешами XML, разрешениями,
spp или AO radius получают отдельные таблицы; PT и AO с одинаковыми сценой,
разрешением и spp объединяются в одной сводке. Повторы одного эксперимента
агрегируются максимумом каждой метрики. Отсутствующая метрика обозначена `--`.
Детальные таблицы по каждому запуску сохраняются ниже сводных.

Измерения RX 6800 XT выполнены 23 ноября 2026 года на RADV Mesa 25.2.8,
`test_035`, 1024×1024, 256 spp; радиус RTAO — 2 единицы сцены.
Для Arc A380 использованы те же сцена, разрешение, spp и радиус, все пять
вариантов и драйвер ANV Mesa 24.3.0-devel (`9fc8668b66`). На этой машине
для захвата выбирался Intel ICD:

```bash
export VK_DRIVER_FILES=/usr/share/vulkan/icd.d/intel_icd.x86_64.json
```

Повторный запуск всех пяти конфигураций Arc A380 воспроизвёл статистики
скомпилированных ядер. Изображения `hw-mega-p` для PT и RTAO темнее, чем у
`hw-mega`, и отличаются между повторными запусками. Для PT сгенерированный
host code делит число проходов на размер подгруппы устройства (32), тогда как
ANV вернул размер 8 для `PathTraceMega`. Данные этого варианта пригодны для
сводки по скомпилированным ядрам. Эквивалентность рендеринга требует проверки.

Сборки, SPIR-V и изображения остаются в игнорируемом `profiling_results`:
на другой машине они создаются командой `build`. Метаданные старых JSON
содержат исходные локальные пути и ревизию HEAD до коммита эксперимента;
`source_dirty: true` указывает на изменения, сохранённые в этой ветке.

`testing/profile_gpu.py` использует Python 3.9+ и стандартную библиотеку.
Скрипт собирает пять независимых вариантов, запускает их и сохраняет исходные
статистики драйвера, параметры эксперимента и LaTeX-таблицы по каждому ядру и
максимумы по ядрам. На каждой машине нужно собирать варианты и запускать
эксперимент локально. Генерация и автоматическая сборка проверены на Linux.

| Вариант | kmake | Трассировка | Исполнение |
|---|---|---|---|
| `sw-wave` | `kmake_wave_swrt_v1.json` | BVH2FatRT | wavefront |
| `sw-mega` | `kmake_mega_swrt_v1.json` | BVH2FatRT | megakernel |
| `hw-wave` | `kmake_wave.json` | Vulkan ray query | wavefront |
| `hw-mega` | `kmake_mega.json` | Vulkan ray query | megakernel |
| `hw-mega-p` | `kmake_mega_persistent.json` | Vulkan ray query | megakernel + persistent rays |

В самой статье также есть SW persistent rays: он не входит в эти пять вариантов.

## Подготовка и сборка

Нужны CMake, C++17, Vulkan headers/loader, `slangc` в `PATH`, собранный
`kernel_slicer`, подмодули Hydra и её обычные зависимости (в том числе Embree/TBB).
SW-ветка использует HeavyRT. В этой старой Hydra он не оформлен как подмодуль.
Проверена следующая ревизия:

```bash
git submodule update --init --recursive
git clone https://github.com/msu-graphics-group/HeavyRT.git external/HeavyRT
git -C external/HeavyRT checkout ee2721c5726aa8841a088f2c85ded2351539e702

python3 testing/profile_gpu.py build \
  --slicer /home/frol/PROG/kernel_slicer
```

`build` автоматически применяет `profiling/heavyrt-compat.patch`, если он ещё
не применён. Патч добавляет `UnderlyingImpl`, аргументы конструктора по умолчанию,
исключает CPU metrics/Timer/allocator из входа слайсера и исправляет пути allocator.
Каталог `core/common` HeavyRT нельзя добавлять в C++ include directories перед
Hydra LiteMath: там находится несовместимая старая копия `LiteMath.h`.

Для существующего клона можно скопировать нужную ревизию в `external/HeavyRT`.
Каталог зависимости игнорируется Git; воспроизводимые изменения сохранены в патче.

Сборки появляются в `profiling_results/builds/<variant>/`: исполняемый файл,
собственные шейдеры, `build.json` с ревизиями и SHA-256, журналы генерации и сборки.
Общий рабочий каталог CMake используется последовательно. Запускать две команды
`build` одновременно в одном клоне нельзя: слайсер перезаписывает generated sources.
Скрипт очищает `shaders_generated` перед каждой генерацией.

Можно выбрать отдельные варианты, число потоков и пересобрать после правок:

```bash
python3 testing/profile_gpu.py build --slicer /path/to/kernel_slicer \
  --variants hw-wave hw-mega hw-mega-p --jobs 8 --rebuild
```

Существующие завершённые сборки пропускаются; после изменения исходников нужен
`--rebuild`. Перенос каталогов со сборками не заменяет сборку на другой машине:
обычные динамические зависимости Hydra должны быть доступны загрузчику.

## Запуск PT и RTAO

Для параметров таблицы 5 по умолчанию используются 1024×1024 и 256 spp.
Пути к сценам задаются явно, чтобы не спутать локальные тестовые сцены со сценами
статьи. Названия Coffee, Spaceship, Staircase и Kitchen — метки в отчёте:

```bash
python3 testing/profile_gpu.py run --gpu-id 0 --integrator mispt \
  --scene Coffee=/path/to/HydraScenes/coffee/statex_00001.xml \
  --scene Spaceship=/path/to/HydraScenes/spaceship/statex_00001.xml \
  --scene Staircase=/path/to/HydraScenes/staircase/statex_00001.xml \
  --scene Kitchen=/path/to/HydraScenes/kitchen/statex_00001.xml \
  --output profiling_results/pt

python3 testing/profile_gpu.py run --gpu-id 0 --integrator rtao --ao-radius 2 \
  --scene Coffee=/path/to/HydraScenes/coffee/statex_00001.xml \
  --output profiling_results/rtao
```

`--scene` повторяется для каждого XML; `--scene-dir` задаёт альтернативный корень
ресурсов сцены. `--variants` ограничивает варианты. `--repeats N` сохраняет каждый
запуск отдельно; для статических регистров обычно достаточно одного запуска.
`--width`, `--height`, `--spp` переопределяют настройки. Для persistent rays число
пикселей и spp должны быть кратны размеру подгруппы. Для RTAO spp должно быть
степенью двойки. Рекомендуются параметры таблицы: 1024×1024, 256 spp.

Быстрая проверка на включённой в репозиторий сцене:

```bash
python3 testing/profile_gpu.py run --scene Smoke=scenes/test_035/statex_00001.xml \
  --width 64 --height 64 --spp 256 --integrator rtao \
  --output profiling_results/smoke-rtao
```

Одна уже собранная Hydra с этим профилировщиком также может быть запущена явно:

```bash
python3 testing/profile_gpu.py run --variants hw-mega-p \
  --binary bin-release/hydra --resource-dir . \
  --scene Smoke=scenes/test_035/statex_00001.xml --output profiling_results/manual
```

При таком запуске соответствие бинарника метке варианта обеспечивает пользователь.
В автоматических сборках проверяются хеши бинарника и шейдеров. Параметр
`-resource_dir` самой Hydra теперь применяется до создания pipeline: каждый вариант
загружает свои SPIR-V, даже когда запускается из корня репозитория.

## Что измеряется

Включается `VK_KHR_pipeline_executable_properties`, запрашивается feature
`pipelineExecutableInfo`, compute pipelines создаются с
`VK_PIPELINE_CREATE_CAPTURE_STATISTICS_BIT_KHR`. После рендеринга считываются
все executable statistics каждого использованного pipeline. Это характеристики
скомпилированного шейдера: не количество переменных GLSL/Slang или SPIR-V IDs.
См. [документацию Khronos](https://docs.vulkan.org/refpages/latest/refpages/source/vkGetPipelineExecutableStatisticsKHR.html).

На RADV доступны SGPR/VGPR, spills, scratch/LDS и `Subgroups per SIMD`.
Последний показатель — статический предел одновременно резидентных подгрупп,
который возвращает драйвер. Он не является achieved occupancy, измеренной во время
исполнения, и не переводится скриптом в проценты. Для achieved occupancy нужны
аппаратные счётчики vendor profiler. Доступные названия и единицы метрик зависят
от GPU и драйвера. В LaTeX сохраняются их исходные названия и определения;
в raw JSON сохраняются все метрики, в том числе не вошедшие в таблицу.

На Arc A380 ANV возвращает ширину SIMD, пиковое регистровое давление
`Max live registers` до аллокации, число операций spill/fill, scratch и SLM.
Скрипт сохраняет фиксированный `LocalSize` из SPIR-V в JSON запуска и вычисляет
теоретический предел по [модели Intel](https://oneapi-src.github.io/oneAPI-samples/Tools/GPU-Occupancy-Calculator/).
Число подгрупп в workgroup округляется вверх по ширине SIMD. В Xe-core помещается
столько целых workgroups, сколько позволяют 128 аппаратных слотов и 128 KiB SLM.
Выделение SLM для workgroup округляется вверх до степени двух с минимумом 1 KiB
при ненулевой потребности. Результат в процентах — доля слотов, занятых на таком
теоретическом пределе. В отчёте для Arc приводится предел архитектуры, а фактическая
загрузка определяется счётчиками во время исполнения. При отсутствии входных
данных для расчёта значение не выводится.

Включаются только compute kernels, для которых записан dispatch во время рендера.
Построение сцены и `PackXY` исключены. `recorded_dispatches` считает запись команд,
а не число фактических исполнений: Hydra может повторно отправлять один command buffer.
`submit_wait_ms` — время submit/wait из Hydra, не аппаратная метрика времени ядра
и не надёжный benchmark под профилировщиком.

Для каждой метрики отдельно выводится максимум по использованным ядрам. Максимумы
VGPR, LDS и residency могут относиться к разным ядрам. Максимум residency не
характеризует occupancy всего wavefront pipeline; для анализа смотрите и строки
тяжёлых ядер. Отсутствующие метрики не подменяются нулём.

Если расширение/feature отсутствует или драйвер не возвращает статистики, скрипт
сохраняет raw JSON и журнал, сообщает ошибку и завершается ненулевым кодом.
Для такой машины нужен поддерживаемый vendor profiler. Чужие GPU и драйверы
должны быть проверены локальным запуском.
Для NVIDIA Vulkan поддерживается
[Shader Profiler в Nsight Graphics](https://docs.nvidia.com/nsight-graphics/UserGuide/shader-profiler.html):
он показывает registers per thread, theoretical warp occupancy и позволяет
экспортировать таблицу в CSV. Это отдельный путь сбора, если стандартное Vulkan
расширение не доступно; данный скрипт такой vendor capture не автоматизирует.

Обычный запуск Hydra не устанавливает Vulkan hooks. Они работают только с
`-profile-json /path/to/profile.json` и рассчитаны на существующий однопоточный
Vulkan host приложения.

## RTAO в Hydra

`-integrator rtao` добавляет в `Integrator` адаптацию HeavyRT
`examples/02_rtao/rtao.cpp`: первичный луч через центр пикселя, интерполяция нормалей,
cosine-weighted AO rays, nearest-hit visibility, 128×128 tile с четырёхмерным
Sobol/Niederreiter bucketing. Таблица создаётся существующим `mlt/rnd_qmc` на CPU.
Результат — линейная видимость [0,1]; промах первичного луча даёт чёрный фон.

В RTAO первичный луч трассируется один раз на пиксель, затем выполняются все AO
samples внутри control function. Wavefront разбивает стадии на отдельные ядра.
Persistent rays распределяют samples пикселя по lanes и выполняют subgroup reduction;
host запускает `AmbientOcclusionBlock(..., 1)`. В PT деление числа запусков на размер
подгруппы уже выполняется generated host code. Дополнительно делить `-spp` не нужно.

Перенесены алгоритм и sampling, но исходники не идентичны старому отдельному
рендеру. Используются Hydra loader, vertex packing, нормали без старого сжатия,
ориентация нормали по первичному лучу и floating-point framebuffer.
`--ao-radius` задаётся **в единицах сцены**, по умолчанию 2. В HeavyRT радиус мог
масштабироваться относительно bounding box: для сравнения с прежними измерениями
установите тот же физический радиус явно. Новые времена не следует представлять
как точное повторение исторических чисел таблицы 5.

SW BVH2FatRT рассчитан на статические triangle meshes. Неподдерживаемая custom
geometry и motion instances завершаются явной ошибкой. HW PT сохраняет исходные
возможности Hydra; RTAO использует статический снимок геометрии.
В SW-сборках включён `REMAP_PRIM_ID`: индексы треугольников после перестановки
BVH переводятся обратно в исходные перед чтением нормалей и материалов Hydra.

Hydra передаёт в pipeline specialization constants, включающие типы материалов
и функции, обнаруженные в сцене. Поэтому регистровые числа, полученные на
`test_035`, характеризуют эту специализацию и не гарантированно совпадают с
Coffee/Spaceship/Staircase/Kitchen. В отчёте сцена указана явно.

## LaTeX и перенос результатов

Каждая команда `run` пишет `profile.tex`, `*.json` с параметрами и raw profile,
`*.raw.json`, журналы и изображения. Для пересоздания общего отчёта можно передать
каталоги запусков, включая скопированные с других GPU:

```bash
python3 testing/profile_gpu.py report \
  profiling_results/pt profiling_results/rtao /path/to/another-gpu/runs \
  --output profiling_results/supplement.tex
```

В документе нужны `\usepackage{longtable}` и `\input{.../supplement.tex}`.
Таблицы рассчитаны на полноширинное приложение/supplement, а не на узкую колонку IEEE.
В отчётах сохраняются модель GPU, vendor/device ID, версии драйвера/API, параметры
рендера, хеш scene XML, ревизии и хеши бинарника/шейдеров. Хеш XML не охватывает
внешние mesh/texture assets; сохраняйте использованные scene libraries отдельно.

Для текста статьи можно использовать формулировку, уточнив значения после запусков:

> Supplementary Material reports per-kernel compiler statistics for the wavefront,
> megakernel, and persistent-rays variants, including register allocation, spills,
> and driver-reported static residency limits where available. These limits are
> distinct from achieved occupancy measured with hardware counters.

Ссылка относится к дополнению с таблицами; утверждение о росте register pressure
следует оставить только если оно подтверждено значениями на соответствующем GPU.
