# Установка и загрузка необходимых пакетов           Dont checked !!!!!
library(sf)
library(terra)
library(stats)
library(FNN)      # Для KNN
library(ggplot2)  # Для визуализации
library(svglite)  # Для сохранения графиков
library(dplyr)    # Для манипуляций с данными

# Очистка памяти
gc()

# Настройки terra для работы с диском
temp_dir <- "D:/temp_terra1"
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
terraOptions(tempdir = temp_dir, memfrac = 0.3)

# Загрузка точек присутствия нового вида (замените путь на ваш файл)
train_balanced <- st_read("new_species_presence.shp")  # Укажите путь к данным нового вида
st_crs(train_balanced) <- "+proj=longlat +datum=WGS84"
cat("CRS train_balanced:", st_crs(train_balanced)$input, "\n")
presence_vect <- vect(train_balanced)

# Загрузка полигона
polygon <- st_read("15_11Layers.shp")  # Если нужен другой полигон, измените путь
if (st_crs(polygon) != st_crs(train_balanced)) {
  polygon <- st_transform(polygon, crs = "+proj=longlat +datum=WGS84")
  cat("CRS полигона приведено:", st_crs(polygon)$input, "\n")
}
polygon_vect <- vect(polygon)

# Загрузка растровых слоев
path_30s <- "t30"  # Если путь отличается, измените
rasters_30s <- list.files(path = path_30s, pattern = "\\.tif$", full.names = TRUE)
if (length(rasters_30s) == 0) {
  stop("Не найдено растровых файлов в папке t30. Проверьте путь.")
}

# Проверка и создание env_layers_all
if (!exists("env_layers_all") || !all(file.exists(sources(env_layers_all)))) {
  cat("env_layers_all недоступен или файл отсутствует. Загружаем заново...\n")
  env_layers_all <- rast(rasters_30s)
  crs(env_layers_all) <- "+proj=longlat +datum=WGS84"
  env_layers_all <- crop(env_layers_all, polygon_vect)
  
  cropped_file_all <- file.path(temp_dir, "cropped_rasters_all.tif")
  if (file.exists(cropped_file_all)) {
    unlink(cropped_file_all, force = TRUE)
  }
  writeRaster(env_layers_all, filename = cropped_file_all, overwrite = TRUE)
  env_layers_all <- rast(cropped_file_all)
  cat("env_layers_all успешно создан и сохранен в", cropped_file_all, "\n")
} else {
  cat("env_layers_all уже доступен.\n")
}

# Выбираем слои для маски SRE (можно изменить список)
selected_layers <- c("wc2.1_30s_bio_11", "wc2.1_30s_bio_15", "wc2.1_30s_bio_3", "wc2.1_30s_bio_9")
raster_names <- gsub(".tif", "", basename(rasters_30s))
selected_indices <- which(raster_names %in% selected_layers)
if (length(selected_indices) != length(selected_layers)) {
  stop("Не удалось найти все указанные слои для маски. Проверьте имена файлов в папке t30.")
}
selected_rasters <- rasters_30s[selected_indices]
env_layers_mask <- rast(selected_rasters)

# Проверка и установка CRS растров для маски
if (is.na(crs(env_layers_mask)) || crs(env_layers_mask) != "+proj=longlat +datum=WGS84") {
  crs(env_layers_mask) <- "+proj=longlat +datum=WGS84"
}
names(env_layers_mask) <- gsub(".tif", "", basename(selected_rasters))
env_layers_mask <- crop(env_layers_mask, polygon_vect)

# Сохранение обрезанных растров для маски
cropped_file_mask <- file.path(temp_dir, "cropped_rasters_mask.tif")
if (file.exists(cropped_file_mask)) {
  unlink(cropped_file_mask, force = TRUE)
}
writeRaster(env_layers_mask, filename = cropped_file_mask, overwrite = TRUE)
env_layers_mask <- rast(cropped_file_mask)

# Извлечение значений для маски SRE
env_presence_mask <- terra::extract(env_layers_mask, presence_vect)[, -1]
if (any(is.na(env_presence_mask))) {
  warning("Найдены NA в данных присутствия для маски. Удаляю строки с NA.")
  valid_rows <- complete.cases(env_presence_mask)
  env_presence_mask <- env_presence_mask[valid_rows, ]
  presence_vect <- presence_vect[valid_rows]
}
gc()

# Вычисление квантилей для SRE
quantiles <- apply(env_presence_mask, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
})

# Создание маски SRE
create_sre_mask <- function(env_layers, quantiles, temp_dir) {
  mask <- env_layers[[1]] * 0 + 1
  for (i in 1:nlyr(env_layers)) {
    layer <- env_layers[[i]]
    q_low <- quantiles[1, i]
    q_high <- quantiles[2, i]
    temp_mask <- layer < q_low | layer > q_high
    temp_file <- file.path(temp_dir, paste0("temp_mask_", i, ".tif"))
    if (file.exists(temp_file)) unlink(temp_file)
    writeRaster(temp_mask, filename = temp_file, overwrite = TRUE)
    mask <- mask * rast(temp_file)
    gc()
  }
  return(mask)
}
sre_mask <- create_sre_mask(env_layers_mask, quantiles, temp_dir)

# Сохранение маски SRE
sre_mask_file <- file.path(temp_dir, "sre_mask.tif")
if (file.exists(sre_mask_file)) unlink(sre_mask_file)
writeRaster(sre_mask, filename = sre_mask_file, overwrite = TRUE)
sre_mask <- rast(sre_mask_file)
sre_mask[sre_mask == 0] <- NA
sre_mask <- mask(sre_mask, polygon_vect)

# Генерация псевдоотсутствий
n_pseudo <- nrow(train_balanced)  # Соотношение 1:1 (можно изменить)
pseudo_absences <- spatSample(sre_mask, size = n_pseudo, method = "random", na.rm = TRUE, as.points = TRUE)

# Извлечение координат и данных для присутствий (все 20 ковариат)
presence_coords <- crds(presence_vect)
env_presence <- terra::extract(env_layers_all, presence_vect)[, -1]
env_presence <- cbind(data.frame(x = presence_coords[, 1], y = presence_coords[, 2]), env_presence)
env_presence$nf <- 1

if (any(is.na(env_presence))) {
  na_count <- sum(!complete.cases(env_presence))
  warning("Обнаружены NA в данных присутствия. Удалено строк: ", na_count)
  env_presence <- na.omit(env_presence)
  presence_vect <- presence_vect[complete.cases(terra::extract(env_layers_all, presence_vect)[, -1]), ]
  presence_coords <- crds(presence_vect)
}

# Извлечение координат и данных для псевдоотсутствий (все 20 ковариат)
pseudo_coords <- crds(pseudo_absences)
env_pseudo <- terra::extract(env_layers_all, pseudo_absences)[, -1]
env_pseudo <- cbind(data.frame(x = pseudo_coords[, 1], y = pseudo_coords[, 2]), env_pseudo)
env_pseudo$nf <- 0

if (any(is.na(env_pseudo))) {
  na_count <- sum(!complete.cases(env_pseudo))
  warning("Обнаружены NA в данных псевдоотсутствий. Удалено строк: ", na_count)
  env_pseudo <- na.omit(env_pseudo)
  pseudo_absences <- pseudo_absences[complete.cases(terra::extract(env_layers_all, pseudo_absences)[, -1]), ]
  pseudo_coords <- crds(pseudo_absences)
}

# Объединение данных
combined_data <- rbind(env_presence, env_pseudo)

# Проверка структуры combined_data
cat("Структура combined_data перед разбиением:\n")
str(combined_data)

# Параметры разбиения (можно настроить)
params <- list(
  test_size = 0.2,
  max_attempts = 90000,
  n_clusters = 12,
  w_cross_dist = 7.0,
  w_coverage = 9.0,
  w_balance = 60.0,
  w_dist_diff = 0.1,
  max_points_per_cluster = 2
)

# Функции для вычисления пространственных метрик
check_spatial_metrics <- function(data) {
  if (nrow(data) < 2) {
    return(list(mean_dist = 0, min_dist = 0, max_dist = 0))
  }
  coords <- as.matrix(data[, c("x", "y")])
  dist_matrix <- as.matrix(dist(coords))
  return(list(
    mean_dist = mean(dist_matrix[upper.tri(dist_matrix)]),
    min_dist = min(dist_matrix[upper.tri(dist_matrix)]),
    max_dist = max(dist_matrix[upper.tri(dist_matrix)])
  ))
}

check_cross_distance <- function(train_data, test_data) {
  if (nrow(train_data) == 0 || nrow(test_data) == 0) {
    return(0)
  }
  knn <- get.knnx(train_data[, c("x", "y")], test_data[, c("x", "y")], k = 1)
  return(mean(knn$nn.dist))
}

# Функция разбиения на основе кластеров
balanced_cluster_split <- function(data, params) {
  if (!all(c("x", "y") %in% colnames(data))) {
    stop("Данные должны содержать колонки 'x' и 'y'")
  }
  if (!is.numeric(data$x) || !is.numeric(data$y)) {
    stop("Столбцы 'x' и 'y' должны быть числовыми")
  }
  coords <- as.matrix(data[, c("x", "y")])
  set.seed(42)
  clusters <- kmeans(coords, centers = params$n_clusters)$cluster
  data$cluster <- as.factor(clusters)
  cat("\nРаспределение всех точек по кластерам:\n")
  print(table(data$cluster))
  target_test <- round(table(data$cluster) * params$test_size)
  target_test[target_test == 0] <- 1
  cluster_counts <- table(data$cluster)
  available_clusters <- names(cluster_counts[cluster_counts > 0])
  best_split <- NULL
  best_score <- -Inf
  for (attempt in 1:params$max_attempts) {
    test_indices <- unlist(lapply(names(target_test), function(cl) {
      cluster_pts <- which(data$cluster == cl)
      sample(cluster_pts, min(target_test[cl], length(cluster_pts)), replace = FALSE)
    }))
    temp_test_data <- data[test_indices, ]
    cluster_coverage <- table(temp_test_data$cluster)
    if (!all(available_clusters %in% names(cluster_coverage))) {
      missing_clusters <- available_clusters[!available_clusters %in% names(cluster_coverage)]
      for (cl in missing_clusters) {
        cluster_pts <- which(data$cluster == cl & !(1:nrow(data) %in% test_indices))
        if (length(cluster_pts) > 0) {
          test_indices <- c(test_indices, sample(cluster_pts, 1))
        }
      }
    }
    temp_test_data <- data[test_indices, ]
    cluster_counts_test <- table(temp_test_data$cluster)
    for (cl in names(cluster_counts_test)) {
      if (cluster_counts_test[cl] > params$max_points_per_cluster) {
        cluster_indices <- which(temp_test_data$cluster == cl)
        excess <- cluster_counts_test[cl] - params$max_points_per_cluster
        remove_indices <- sample(cluster_indices, excess)
        test_indices <- test_indices[!test_indices %in% temp_test_data$index[remove_indices]]
      }
    }
    data$index <- 1:nrow(data)
    temp_test_data <- data[test_indices, ]
    cluster_coverage <- table(temp_test_data$cluster)
    if (!all(available_clusters %in% names(cluster_coverage))) {
      missing_clusters <- available_clusters[!available_clusters %in% names(cluster_coverage)]
      for (cl in missing_clusters) {
        cluster_pts <- which(data$cluster == cl & !(1:nrow(data) %in% test_indices))
        if (length(cluster_pts) > 0) {
          test_indices <- c(test_indices, sample(cluster_pts, 1))
        }
      }
    }
    train_data <- data[-test_indices, ]
    test_data <- data[test_indices, ]
    train_data$index <- NULL
    test_data$index <- NULL
    data$index <- NULL
    train_metrics <- check_spatial_metrics(train_data)
    test_metrics <- check_spatial_metrics(test_data)
    cross_dist <- check_cross_distance(train_data, test_data)
    coverage <- mean(table(test_data$cluster) > 0)
    balance <- 1 - sd(table(test_data$cluster)) / mean(table(test_data$cluster))
    dist_diff <- abs(train_metrics$mean_dist - test_metrics$mean_dist)
    score <- (cross_dist * params$w_cross_dist + 
                coverage * params$w_coverage + 
                balance * params$w_balance - 
                dist_diff * params$w_dist_diff)
    if (score > best_score) {
      best_split <- list(train = train_data, test = test_data)
      best_score <- score
    }
  }
  return(best_split)
}

# Применение функции разбиения
set.seed(42)
split_result <- balanced_cluster_split(combined_data, params)

# Разделение на тренировочные и тестовые наборы
train_data <- split_result$train
test_data <- split_result$test

# Вычисление метрик
train_metrics <- check_spatial_metrics(train_data)
test_metrics <- check_spatial_metrics(test_data)
cross_dist <- check_cross_distance(train_data, test_data)
coverage <- length(unique(test_data$cluster)) / params$n_clusters

# Подсчет доли присутствий и псевдоотсутствий в наборах
train_presence_ratio <- mean(train_data$nf == 1) * 100
test_presence_ratio <- mean(test_data$nf == 1) * 100

# Визуализация
knn <- get.knnx(train_data[, c("x", "y")], test_data[, c("x", "y")], k = 1)
segment_data <- data.frame(
  x = test_data$x,
  y = test_data$y,
  xend = train_data$x[knn$nn.index],
  yend = train_data$y[knn$nn.index]
)

final_plot <- ggplot() +
  geom_point(data = combined_data, aes(x, y, color = as.factor(nf)), size = 2, alpha = 0.3) +
  geom_point(data = train_data, aes(x, y, color = "Training"), size = 3, alpha = 0.7) +
  geom_point(data = test_data, aes(x, y, color = "Testing"), shape = 17, size = 4) +
  geom_segment(
    data = segment_data,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "darkgrey", linetype = "dashed", alpha = 0.5
  ) +
  scale_color_manual(
    name = "Dataset",
    values = c("0" = "grey90", "1" = "grey50", "Training" = "#3575b5", "Testing" = "#e74c3c"),
    labels = c("0" = "Pseudo-absence", "1" = "Presence", "Training" = "Training", "Testing" = "Testing")
  ) +
  labs(
    title = "Balanced Cluster-Based Data Split",
    subtitle = sprintf(
      "Training: %d points (%.1f%% presence) | Testing: %d points (%.1f%% presence) | Cross-dist: %.2f | Cluster coverage: %.0f%%",
      nrow(train_data), train_presence_ratio, nrow(test_data), test_presence_ratio, cross_dist, coverage * 100
    ),
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Сохранение графика
ggsave("balanced_cluster_split_plot.png", final_plot, width = 8, height = 6, dpi = 300)
ggsave("balanced_cluster_split_plot.tiff", final_plot, width = 8, height = 6, dpi = 300, compression = "lzw")
ggsave("balanced_cluster_split_plot.pdf", final_plot, width = 8, height = 6, dpi = 300)

# Сохранение данных
write.table(test_data, file = "test_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(train_data, file = "train_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Создание текстового отчета
report_text <- paste(
  "=== SPLITTING RESULTS ===",
  "\nParameters:",
  sprintf("- Test size: %.1f", params$test_size),
  sprintf("- Max attempts: %d", params$max_attempts),
  sprintf("- Number of clusters: %d", params$n_clusters),
  sprintf("- Cross-distance weight: %.1f", params$w_cross_dist),
  sprintf("- Coverage weight: %.1f", params$w_coverage),
  sprintf("- Balance weight: %.1f", params$w_balance),
  sprintf("- Distance difference weight: %.1f", params$w_dist_diff),
  sprintf("- Max points per cluster: %d", params$max_points_per_cluster),
  "\nDataset sizes:",
  sprintf("- Training set: %d points (%.1f%%)", nrow(train_data), nrow(train_data)/nrow(combined_data)*100),
  sprintf("- Testing set: %d points (%.1f%%)", nrow(test_data), nrow(test_data)/nrow(combined_data)*100),
  "\nPresence/Absence distribution:",
  sprintf("- Training set presence: %.1f%%", train_presence_ratio),
  sprintf("- Testing set presence: %.1f%%", test_presence_ratio),
  "\nSpatial metrics:",
  sprintf("- Mean distance in train: %.3f", train_metrics$mean_dist),
  sprintf("- Mean distance in test: %.3f", test_metrics$mean_dist),
  sprintf("- Cross-group distance: %.3f", cross_dist),
  sprintf("- Cluster coverage: %.1f%%", coverage * 100),
  "\nTesting points per cluster:",
  capture.output(print(table(test_data$cluster))),
  sep = "\n"
)

# Сохранение отчета
writeLines(report_text, "splitting_report.txt")

# Вывод графика и отчета
print(final_plot)
cat(report_text)