#!/usr/bin/env python3
import sys
import os
import csv
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import concurrent.futures

def compute_percentage(record):
    """
    Вычисляет процент CDS, где значение 'gene' совпадает с 'locus_tag'.
    Возвращает кортеж (процент совпадения, процент несовпадения).
    """
    total = 0
    match = 0
    for feature in record.features:
        if feature.type != "CDS":
            continue
        lt = feature.qualifiers.get("locus_tag", [None])[0]
        gene = feature.qualifiers.get("gene", [lt])[0]
        if lt is not None:
            total += 1
            if gene == lt:
                match += 1
    return (match / total * 100, 100 - match / total * 100) if total > 0 else (0, 0)

def compute_identity(target_seq, ref_seq):
    """
    Вычисляет идентичность между двумя белковыми последовательностями с использованием
    глобального выравнивания через Align.PairwiseAligner.
    
    В данном случае используется метод score(), который возвращает оптимальный счёт выравнивания.
    Идентичность вычисляется как отношение best_score / max(len(target_seq), len(ref_seq)).
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    best_score = aligner.score(target_seq, ref_seq)
    max_len = max(len(target_seq), len(ref_seq))
    return best_score / max_len if max_len > 0 else 0

def find_best_match(target_translation, ref_features, threshold=0.7):
    """
    Ищет среди CDS референсного файла наиболее подходящую запись по выравниванию
    белковых последовательностей (translation). Возвращает (ref_feature, identity_ratio)
    или (None, 0), если подходящего совпадения нет.
    """
    best_ratio = 0
    best_feature = None
    for ref in ref_features:
        ref_translation = ref.get("translation")
        if not ref_translation or not target_translation:
            continue
        ratio = compute_identity(target_translation, ref_translation)
        if ratio >= threshold and ratio > best_ratio:
            best_ratio = ratio
            best_feature = ref
    return best_feature, best_ratio

def process_feature(feature_dict, ref_features, threshold=0.7):
    """
    Обрабатывает одну CDS из целевого файла.
    
    Если белковые последовательности (translation) совпадают с идентичностью ≥ threshold,
    то:
      - Если поле gene отсутствует или равно locus_tag, обновляем gene из референсного CDS.
      - Если в поле product целевого CDS не содержится уже обновлённое имя gene,
        то переписываем product из референсного CDS.
        
    Возвращает словарь с обновлёнными значениями и флагами изменений.
    """
    original_gene = feature_dict["gene"]
    locus = feature_dict["locus"]
    original_product = feature_dict["product"]
    target_translation = feature_dict["translation"]

    updated_gene = original_gene
    updated_product = original_product
    change_gene = False
    change_product = False

    if not target_translation:
        return {"global_index": feature_dict["global_index"],
                "updated_gene": updated_gene,
                "updated_product": updated_product,
                "change_gene": change_gene,
                "change_product": change_product,
                "original_gene": original_gene,
                "original_product": original_product}

    best_ref, ratio = find_best_match(target_translation, ref_features, threshold)
    if best_ref is None:
        return {"global_index": feature_dict["global_index"],
                "updated_gene": updated_gene,
                "updated_product": updated_product,
                "change_gene": change_gene,
                "change_product": change_product,
                "original_gene": original_gene,
                "original_product": original_product}

    # Условие (1): если gene отсутствует или равно locus_tag, обновляем gene из референсного CDS
    if (original_gene is None or original_gene == locus):
        ref_gene = best_ref.get("gene")
        ref_locus = best_ref.get("locus")
        if ref_gene is not None and ref_gene != ref_locus:
            updated_gene = ref_gene
            change_gene = True

    # Корректировка product: если в целевом product не содержится уже обновлённое имя gene,
    # то переписываем product из референсного CDS, если оно задано.
    product_lower = original_product.lower() if original_product else ""
    gene_in_product = updated_gene and updated_gene.lower() in product_lower
    if not gene_in_product and best_ref.get("product"):
        new_product = best_ref.get("product")
        if new_product != original_product:
            updated_product = new_product
            change_product = True

    return {"global_index": feature_dict["global_index"],
            "updated_gene": updated_gene,
            "updated_product": updated_product,
            "change_gene": change_gene,
            "change_product": change_product,
            "original_gene": original_gene,
            "original_product": original_product}

def main():
    if len(sys.argv) != 3:
        print("Использование: python скрипт.py целевой.gb референсный.gb")
        sys.exit(1)

    target_file, ref_file = sys.argv[1], sys.argv[2]

    try:
        target_record = SeqIO.read(target_file, "genbank")
        ref_record = SeqIO.read(ref_file, "genbank")
    except Exception as e:
        sys.exit(f"Ошибка чтения файлов: {e}")

    target_match, target_non = compute_percentage(target_record)
    ref_match, ref_non = compute_percentage(ref_record)
    print(f"Целевой файл до обработки: {target_match:.2f}% совпадение gene == locus_tag, {target_non:.2f}% не совпадает")
    print(f"Референсный файл: {ref_match:.2f}% совпадение gene == locus_tag, {ref_non:.2f}% не совпадает")

    # Извлечение CDS из референсного файла
    ref_features = []
    for feat in ref_record.features:
        if feat.type == "CDS":
            ref_features.append({
                "locus": feat.qualifiers.get("locus_tag", [None])[0],
                "gene": feat.qualifiers.get("gene", [None])[0],
                "product": feat.qualifiers.get("product", [None])[0],
                "translation": feat.qualifiers.get("translation", [None])[0]
            })

    # Извлечение CDS из целевого файла и сохранение их глобального индекса
    target_features = []
    for i, feat in enumerate(target_record.features):
        if feat.type != "CDS":
            continue
        target_features.append({
            "global_index": i,  # глобальный индекс в target_record.features
            "locus": feat.qualifiers.get("locus_tag", [None])[0],
            "gene": feat.qualifiers.get("gene", [feat.qualifiers.get("locus_tag", [None])[0]])[0],
            "product": feat.qualifiers.get("product", [None])[0],
            "translation": feat.qualifiers.get("translation", [None])[0]
        })

    # Параллельная обработка CDS; результаты сохраняются в словарь по их глобальному индексу
    results = {}
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_feature, feat, ref_features, 0.7): feat["global_index"] for feat in target_features}
        for future in concurrent.futures.as_completed(futures):
            res = future.result()
            results[res["global_index"]] = res

    # Обновление аннотации и формирование лога изменений, используя глобальные индексы
    log_rows = []
    for global_index, feat_update in results.items():
        feature = target_record.features[global_index]
        if feat_update["updated_gene"]:
            feature.qualifiers["gene"] = [feat_update["updated_gene"]]
        if feat_update["updated_product"]:
            feature.qualifiers["product"] = [feat_update["updated_product"]]
        orig_str = f"gene: {feat_update['original_gene']} | product: {feat_update['original_product']}"
        changes = []
        if feat_update["change_gene"]:
            changes.append(f"gene: {feat_update['original_gene']} -> {feat_update['updated_gene']}")
        if feat_update["change_product"]:
            changes.append(f"product: {feat_update['original_product']} -> {feat_update['updated_product']}")
        change_str = "; ".join(changes)
        log_rows.append({
            "Целевой до изменений": orig_str,
            "Изменения из референсного": change_str,
            "Изменено название гена": "Да" if feat_update["change_gene"] else "Нет",
            "Изменено название продукта": "Да" if feat_update["change_product"] else "Нет"
        })

    updated_match, updated_non = compute_percentage(target_record)
    print(f"Целевой файл после обработки: {updated_match:.2f}% совпадение gene == locus_tag, {updated_non:.2f}% не совпадает")

    # Сохранение обновлённого GenBank‑файла
    out_filename = "edited.gb"
    with open(out_filename, "w") as outfile:
        SeqIO.write(target_record, outfile, "genbank")
    print(f"Отредактированный файл сохранён: {out_filename}")

    # Сохранение лога в CSV
    csv_filename = "results.csv"
    with open(csv_filename, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = ["Целевой до изменений", "Изменения из референсного",
                      "Изменено название гена", "Изменено название продукта"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in log_rows:
            writer.writerow(row)
    print(f"Лог CSV сохранён: {csv_filename}")

    # Если установлен pandas – сохраняем также Excel-лог
    try:
        import pandas as pd
        df = pd.DataFrame(log_rows)
        xlsx_filename = "results.xlsx"
        df.to_excel(xlsx_filename, index=False)
        print(f"Лог Excel сохранён: {xlsx_filename}")
    except ImportError:
        print("pandas не найден, Excel-лог (results.xlsx) не создан.")

if __name__ == '__main__':
    main()
