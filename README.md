# Compgen-2025: A Journey into Computational Genomics

This repository is structured into **two core parts**:

1. 🔬 **Hands-on Projects** — This section showcases several practical analyses where I applied the knowledge from the course to solve real-world problems using modern genomics datasets. These projects reflect my ability to work with high-dimensional data, build analytical workflows, and extract biological insights across diverse contexts such as spatial transcriptomics, multi-omics integration, and predictive modeling.

2. 📚 **Course Notes & Supplementary Materials** — I've also compiled a set of comprehensive notes covering the theoretical foundations taught throughout the course. These notes emphasize core concepts, reasoning behind analytical steps, and relevant background knowledge. Where helpful, I’ve included carefully curated supplementary materials (e.g., papers, reviews, tutorials)—**every one of them personally reviewed** to ensure relevance and clarity.

Whether you're here to explore full-scale genomics projects or to deepen your understanding of computational methods, I hope this repo serves as both a learning resource and a showcase of practical capabilities.

**Course Notes URL:** [https://dewey-wang.github.io/Compgen-2025/](https://dewey-wang.github.io/Compgen-2025/)

---

#  **Hands-on Projects**

## Visium Analysis

-   The markdown is using **10X Visium data (Mouse Brain)** and covers normalization, dimensionality reduction, clustering, visualization, and integration with single-cell data for deconvolution and cell type inference.
-   [Code](https://dewey-wang.github.io/Compgen-2025/module2/Visium-analysis.html)
  
## Xenium Analysis

-   Xenium is a **cell-level spatial transcriptomics platform** that offers subcellular resolution. This allows us to examine not just spatial zones of gene expression but **individual cells and their microenvironments**. Here, we perform a full analysis including: data filtering, disk-based representation for large-scale handling, clustering, marker analysis, and hotspot detection.
-   [Code](https://dewey-wang.github.io/Compgen-2025/module2/Xenium-analysis.html)
  
## Xenium Analysis with Image Registration

-   Registered multimodal images (cell assays, molecular assays, H&E, ROI) to integrate spatial and morphological features. Studied spatially resolved gene expression and microenvironment interactions.
-   [Code](https://dewey-wang.github.io/Compgen-2025/module2/Image-Registration.html)
  
## Breast Cancer Subtype Prediction

-   Integrated METABRIC multi-omics data (gex, cna, clinical); selected features with Laplacian score; trained model with Bayesian hyperparameter tuning; evaluated classification performance.
-   [Code](code/module%203/hw/day%201/noteboook/day1_hw_brca_subtypes.ipynb)

## Survival Analysis

-   Integrated TCGA LGG/GBM data (mutation, cna, clinical); performed EDA; built survival models targeting OS STATUS/OS MONTHS; computed top 5 features with Integrated Gradients; built Cox-PH model with markers and clinical variables; visualized log hazard ratios.
-   [Code](code/module%203/hw/day%202/noteboook/day2_hw_lgg_gbm.ipynb)

## Drug Response Prediction

-   Integrated CCLE/GDSC data (mutation, rna) to build models predicting Erlotinib responses; enhanced neural network models by incorporating STRINGDB protein network with Bayesian hyperparameter tuning; computed top 10 features with Integrated Gradients for correlation analysis.
-   [Code](code/module%203/hw/day%203/noteboook/day3.ipynb)

---


## Course Notes & Supplementary Materials

The theoretical side of this course is just as rich. In this section, I’ve organized detailed notes covering the key ideas from each module, with a strong focus on *why* each method is used—not just *how*. These notes are meant to clarify core statistical and computational genomics concepts, which are crucial for understanding and implementing advanced analyses.

For each module, I also include **supplementary materials**—carefully selected papers, reviews, and tools that extend the core topics. Every supplementary item was **personally reviewed and curated** by me to ensure it truly adds value to the learning process.

#### Module 1 Topics
- Statistics in Genomics  
- Unsupervised Learning  
- Supervised Learning  

#### Module 2 Topics
- Introduction to Spatial Transcriptomics  
- Spatial Transcriptomics Analysis  
- Multiomics Integration  

#### Module 3 Topics
- Multiomics Integration & Neural Network Models  
- Neural Network Training  
- Genomics Model Architectures  


---

## How to Use This Repository

- **Notes**: Browse the HTML files or visit the hosted site ([https://dewey-wang.github.io/Compgen-2025/](https://dewey-wang.github.io/Compgen-2025/)) for a polished reading experience.
- **Code**: Check the code directories for scripts where I apply the knowledge from the course to practical problems in genomics. These include implementations of statistical methods, machine learning models, and spatial transcriptomics workflows.
- **Contributions**: Feel free to suggest improvements, share feedback, or collaborate by opening an issue or pull request!

---

## Acknowledgments

A heartfelt thank you to **Dr. Altuna Akalın**, **Dr. Artür Manukyan**, **Dr. Bora Uyar**, and the entire team at the **Max Delbrück Center for Molecular Medicine (MDC)** for their exceptional teaching and support throughout the course.

---

## Connect With Me

For any questions, issues, or thoughts, please contact via Issues or my email:  
**Ding Yang Wang**  
Email: [deweywang2000@gmail.com](deweywang2000@gmail.com)

