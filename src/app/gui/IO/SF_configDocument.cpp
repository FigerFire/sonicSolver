/// @file SF_configDocument.cpp
/// @brief GUI 工程与中立配置文档的读写实现。

#include "SF_configDocument.h"

#include <QFile>
#include <QFileInfo>
#include <QHash>
#include <QRegularExpression>

#include <algorithm>

namespace SF::GUI {
namespace {

QString leadingWhitespace(const QString& text)
{
    int count = 0;
    while (count < text.size() && text.at(count).isSpace()) ++count;
    return text.left(count);
}

bool normalizeCoordinateTuple(const QString& text, QString* normalized)
{
    const QString number = QStringLiteral(
        R"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)");
    const QRegularExpression pattern(
        QStringLiteral(R"(^\s*\(?\s*(%1)\s*[,\s]+(%1)\s*[,\s]+(%1)\s*\)?\s*$)")
            .arg(number));
    const QRegularExpressionMatch match = pattern.match(text);
    if (!match.hasMatch()) return false;
    if (normalized) {
        *normalized = QStringLiteral("(%1 %2 %3)")
                          .arg(match.captured(1), match.captured(2),
                               match.captured(3));
    }
    return true;
}

QString stripFoamComments(QString text)
{
    text.remove(QRegularExpression(QStringLiteral(R"(/\*[\s\S]*?\*/)")));
    text.remove(QRegularExpression(QStringLiteral(R"(//[^\n]*)")));
    return text;
}

int findMatching(const QStringList& lines,
                 int startLine,
                 QChar open,
                 QChar close)
{
    int depth = 0;
    for (int i = startLine; i < lines.size(); ++i) {
        const QString clean = stripFoamComments(lines.at(i));
        for (const QChar ch : clean) {
            if (ch == open) ++depth;
            else if (ch == close) {
                --depth;
                if (depth == 0) return i;
            }
        }
    }
    return -1;
}

void addLineEntry(QVector<ConfigEntry>* entries,
                  const QString& section,
                  const QString& key,
                  const QString& value,
                  int line)
{
    ConfigEntry entry;
    entry.section = section;
    entry.sectionLabel = section;
    entry.key = key;
    entry.value = value;
    entry.startLine = line;
    entry.endLine = line;
    entry.caseSection = true;
    entries->push_back(entry);
}

bool isBlockMeshDict(const QString& path)
{
    return QFileInfo(path).fileName().compare(
               QStringLiteral("blockMeshDict"), Qt::CaseInsensitive) == 0;
}

int findBlockMeshListEnd(const QStringList& lines, const QString& listName)
{
    const QRegularExpression listPattern(
        QStringLiteral(R"(^\s*%1\s*(?:\(|$))")
            .arg(QRegularExpression::escape(listName)),
        QRegularExpression::CaseInsensitiveOption);
    for (int i = 0; i < lines.size(); ++i) {
        if (!listPattern.match(stripFoamComments(lines.at(i))).hasMatch()) continue;
        return findMatching(lines, i, QLatin1Char('('), QLatin1Char(')'));
    }
    return -1;
}

bool patchRange(const QStringList& lines, const QString& patchName,
                int* openLine, int* closeLine)
{
    const int boundaryEnd = findBlockMeshListEnd(lines, QStringLiteral("boundary"));
    if (boundaryEnd < 0) return false;
    const QRegularExpression namePattern(
        QStringLiteral(R"(^\s*%1\s*(?:\{|$))")
            .arg(QRegularExpression::escape(patchName)));
    for (int line = 0; line < boundaryEnd; ++line) {
        if (!namePattern.match(stripFoamComments(lines.at(line))).hasMatch()) continue;
        int open = line;
        if (!stripFoamComments(lines.at(open)).contains(QLatin1Char('{'))) ++open;
        if (open >= boundaryEnd
            || !stripFoamComments(lines.at(open)).contains(QLatin1Char('{'))) continue;
        const int close = findMatching(lines, open, QLatin1Char('{'), QLatin1Char('}'));
        if (close < 0) return false;
        if (openLine) *openLine = open;
        if (closeLine) *closeLine = close;
        return true;
    }
    return false;
}

} // namespace

bool ConfigEntry::isBoolean() const
{
    const QString normalized = value.trimmed().toLower();
    return normalized == QStringLiteral("true")
        || normalized == QStringLiteral("false");
}

bool ConfigDocument::load(const QString& path, QString* error)
{
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (error) *error = file.errorString();
        return false;
    }

    path_ = QFileInfo(path).absoluteFilePath();
    QString content = QString::fromUtf8(file.readAll());
    content.replace(QStringLiteral("\r\n"), QStringLiteral("\n"));
    content.replace(QLatin1Char('\r'), QLatin1Char('\n'));
    lines_ = content.split(QLatin1Char('\n'), Qt::KeepEmptyParts);
    if (!lines_.isEmpty() && lines_.last().isEmpty()) lines_.removeLast();
    dirty_ = false;
    parse();
    return true;
}

bool ConfigDocument::save(QString* error) const
{
    QFile file(path_);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) {
        if (error) *error = file.errorString();
        return false;
    }

    QByteArray data = lines_.join(QLatin1Char('\n')).toUtf8();
    data.append('\n');
    if (file.write(data) != data.size()) {
        if (error) *error = file.errorString();
        return false;
    }
    return true;
}

bool ConfigDocument::setValue(int entryIndex,
                              const QString& value,
                              QString* error)
{
    if (entryIndex < 0 || entryIndex >= entries_.size()) {
        if (error) *error = QStringLiteral("无效的配置项索引");
        return false;
    }

    const ConfigEntry entry = entries_.at(entryIndex);
    if (entry.startLine < 0 || entry.endLine >= lines_.size()) {
        if (error) *error = QStringLiteral("配置项行范围无效");
        return false;
    }

    QStringList replacement;
    if (entry.coordinateTuple) {
        QString coordinate;
        if (!normalizeCoordinateTuple(value, &coordinate)) {
            if (error) {
                *error = QStringLiteral(
                    "顶点坐标必须是三个数，例如：0 1 2");
            }
            return false;
        }
        replacement << leadingWhitespace(lines_.at(entry.startLine))
                    + coordinate;
    } else if (entry.foamAssignment) {
        QString foamValue = value.trimmed();
        while (foamValue.endsWith(QLatin1Char(';'))) foamValue.chop(1);
        if (foamValue.isEmpty()) {
            if (error) *error = QStringLiteral("OpenFOAM 字典值不能为空");
            return false;
        }
        const QString original = lines_.at(entry.startLine);
        const QString indent = leadingWhitespace(original);
        const QString comment = inlineComment(original);
        QString line = indent + entry.key + QLatin1Char(' ') + foamValue
            + QLatin1Char(';');
        if (!comment.isEmpty()) line += QLatin1Char(' ') + comment;
        replacement << line;
    } else if (entry.caseSection) {
        replacement = value.split(QLatin1Char('\n'), Qt::KeepEmptyParts);
        const QString indent = leadingWhitespace(lines_.at(entry.startLine));
        for (QString& line : replacement) {
            if (!line.isEmpty() && leadingWhitespace(line).isEmpty()) {
                line.prepend(indent);
            }
        }
    } else {
        const QString original = lines_.at(entry.startLine);
        const QString indent = leadingWhitespace(original);
        const QString comment = inlineComment(original);
        QStringList valueLines = value.split(QLatin1Char('\n'), Qt::KeepEmptyParts);
        if (valueLines.isEmpty()) valueLines << QString();

        QString first = indent + entry.key + QStringLiteral(" = ")
            + valueLines.takeFirst().trimmed();
        if (!comment.isEmpty()) first += QLatin1Char(' ') + comment;
        replacement << first;
        for (const QString& line : valueLines) replacement << line;
    }

    for (int line = entry.endLine; line >= entry.startLine; --line) {
        lines_.removeAt(line);
    }
    for (int i = 0; i < replacement.size(); ++i) {
        lines_.insert(entry.startLine + i, replacement.at(i));
    }

    dirty_ = true;
    parse();
    return true;
}

bool ConfigDocument::appendAssignments(const QString& section,
                                       const QStringList& assignments,
                                       bool newSection,
                                       QString* error)
{
    if (section.trimmed().isEmpty() || assignments.isEmpty()) {
        if (error) *error = QStringLiteral("section 或新增内容为空");
        return false;
    }

    if (QFileInfo(path_).fileName().compare(
            QStringLiteral("blockMeshDict"), Qt::CaseInsensitive) == 0) {
        const QString listName = section.trimmed();
        const QRegularExpression listPattern(
            QStringLiteral(R"(^\s*%1\s*(?:\(|$))")
                .arg(QRegularExpression::escape(listName)),
            QRegularExpression::CaseInsensitiveOption);
        int listStart = -1;
        int openLine = -1;
        int depth = 0;
        for (int i = 0; i < lines_.size(); ++i) {
            const QString clean = stripFoamComments(lines_.at(i));
            if (listStart < 0) {
                if (!listPattern.match(clean).hasMatch()) continue;
                listStart = i;
            }
            for (const QChar ch : clean) {
                if (ch == QLatin1Char('(')) {
                    if (openLine < 0) openLine = i;
                    ++depth;
                } else if (ch == QLatin1Char(')') && openLine >= 0) {
                    --depth;
                    if (depth == 0) {
                        QStringList insertion;
                        for (const QString& assignment : assignments) {
                            const int equal = assignment.indexOf(QLatin1Char('='));
                            QString value = equal >= 0
                                ? assignment.mid(equal + 1).trimmed()
                                : assignment.trimmed();
                            if (listName.compare(QStringLiteral("vertices"),
                                                 Qt::CaseInsensitive) == 0) {
                                QString coordinate;
                                if (!normalizeCoordinateTuple(value, &coordinate)) {
                                    if (error) {
                                        *error = QStringLiteral(
                                            "顶点坐标必须是三个数，例如：0 1 2");
                                    }
                                    return false;
                                }
                                value = coordinate;
                            }
                            if (!value.isEmpty()) {
                                insertion << QStringLiteral("    ") + value;
                            }
                        }
                        if (insertion.isEmpty()) {
                            if (error) *error = QStringLiteral("没有可写入的 blockMesh 条目");
                            return false;
                        }
                        for (int offset = 0; offset < insertion.size(); ++offset) {
                            lines_.insert(i + offset, insertion.at(offset));
                        }
                        dirty_ = true;
                        parse();
                        return true;
                    }
                }
            }
        }
        if (error) {
            *error = QStringLiteral("blockMeshDict 中未找到 %1 列表")
                .arg(listName);
        }
        return false;
    }

    QStringList formatted;
    for (QString assignment : assignments) {
        assignment = assignment.trimmed();
        if (!assignment.isEmpty()) {
            formatted.push_back(QStringLiteral("    ") + assignment);
        }
    }
    if (formatted.isEmpty()) {
        if (error) *error = QStringLiteral("没有可写入的配置项");
        return false;
    }

    if (newSection) {
        if (!lines_.isEmpty() && !lines_.last().trimmed().isEmpty()) {
            lines_.push_back(QString());
        }
        lines_.push_back(QStringLiteral("[%1]").arg(section));
        lines_.append(formatted);
    } else {
        const QRegularExpression sectionPattern(
            QStringLiteral(R"(^\s*\[([^\]]+)\]\s*(?:#.*)?$)"));
        int insertion = -1;
        bool inSection = false;
        for (int i = 0; i < lines_.size(); ++i) {
            const QRegularExpressionMatch match =
                sectionPattern.match(lines_.at(i));
            if (match.hasMatch()) {
                if (inSection) break;
                inSection = match.captured(1).trimmed().compare(
                    section, Qt::CaseInsensitive) == 0;
                if (inSection) insertion = i + 1;
                continue;
            }
            if (inSection) insertion = i + 1;
        }
        if (insertion < 0) {
            if (!lines_.isEmpty() && !lines_.last().trimmed().isEmpty()) {
                lines_.push_back(QString());
            }
            lines_.push_back(QStringLiteral("[%1]").arg(section));
            insertion = lines_.size();
        } else if (insertion > 0
                   && !lines_.at(insertion - 1).trimmed().isEmpty()) {
            formatted.prepend(QString());
        }
        for (int i = 0; i < formatted.size(); ++i) {
            lines_.insert(insertion + i, formatted.at(i));
        }
    }

    dirty_ = true;
    parse();
    return true;
}

bool ConfigDocument::addBlockMeshPatch(const QString& name,
                                       const QString& type,
                                       QString* error)
{
    if (!isBlockMeshDict(path_)) {
        if (error) *error = QStringLiteral("只有 blockMeshDict 支持 patch 编辑");
        return false;
    }
    const QString normalizedName = name.trimmed();
    if (!QRegularExpression(QStringLiteral(R"(^[A-Za-z_][A-Za-z0-9_.-]*$)"))
             .match(normalizedName).hasMatch()) {
        if (error) *error = QStringLiteral("patch 名称只能包含字母、数字、_、. 和 -，且不能以数字开头");
        return false;
    }
    int unusedOpen = -1;
    if (patchRange(lines_, normalizedName, &unusedOpen, nullptr)) {
        if (error) *error = QStringLiteral("patch %1 已存在").arg(normalizedName);
        return false;
    }
    const QString normalizedType = type.trimmed().toLower();
    if (normalizedType != QStringLiteral("patch")
        && normalizedType != QStringLiteral("wall")) {
        if (error) *error = QStringLiteral("patch 类型只能是 patch 或 wall");
        return false;
    }
    const int end = findBlockMeshListEnd(lines_, QStringLiteral("boundary"));
    if (end < 0) {
        if (error) *error = QStringLiteral("blockMeshDict 中未找到 boundary 列表");
        return false;
    }
    const QStringList text{
        QStringLiteral("    %1").arg(normalizedName),
        QStringLiteral("    {"),
        QStringLiteral("        type %1;").arg(normalizedType),
        QStringLiteral("        faces ();"),
        QStringLiteral("    }")};
    for (int i = 0; i < text.size(); ++i) lines_.insert(end + i, text.at(i));
    dirty_ = true;
    parse();
    return true;
}

bool ConfigDocument::addBlockMeshPatchFace(const QString& patchName,
                                           const QString& face,
                                           QString* error)
{
    if (!isBlockMeshDict(path_)) {
        if (error) *error = QStringLiteral("只有 blockMeshDict 支持 patch face 编辑");
        return false;
    }
    const QRegularExpression facePattern(
        QStringLiteral(R"(^\s*\(\s*\d+\s+\d+\s+\d+\s+\d+\s*\)\s*$)"));
    if (!facePattern.match(face).hasMatch()) {
        if (error) *error = QStringLiteral("patch 面必须由四个顶点索引组成");
        return false;
    }
    int open = -1;
    int close = -1;
    if (!patchRange(lines_, patchName, &open, &close)) {
        if (error) *error = QStringLiteral("未找到 patch %1").arg(patchName);
        return false;
    }
    const QRegularExpression facesPattern(QStringLiteral(R"(^\s*faces\s*(?:\(|$))"));
    for (int line = open + 1; line < close; ++line) {
        if (!facesPattern.match(stripFoamComments(lines_.at(line))).hasMatch()) continue;
        if (stripFoamComments(lines_.at(line)).contains(
                QRegularExpression(QStringLiteral(R"(faces\s*\(\s*\)\s*;)")))) {
            const QString indent = leadingWhitespace(lines_.at(line));
            lines_.removeAt(line);
            lines_.insert(line, indent + QStringLiteral("faces"));
            lines_.insert(line + 1, indent + QStringLiteral("("));
            lines_.insert(line + 2, indent + QStringLiteral("    %1").arg(face.trimmed()));
            lines_.insert(line + 3, indent + QStringLiteral(");"));
            dirty_ = true;
            parse();
            return true;
        }
        const int facesEnd = findMatching(lines_, line, QLatin1Char('('), QLatin1Char(')'));
        if (facesEnd < 0 || facesEnd > close) break;
        lines_.insert(facesEnd, QStringLiteral("            %1").arg(face.trimmed()));
        dirty_ = true;
        parse();
        return true;
    }
    if (error) *error = QStringLiteral("patch %1 中未找到 faces 列表").arg(patchName);
    return false;
}

bool ConfigDocument::renameBlockMeshPatch(const QString& oldName,
                                          const QString& newName,
                                          QString* error)
{
    if (!isBlockMeshDict(path_)) {
        if (error) *error = QStringLiteral("只有 blockMeshDict 支持 patch 编辑");
        return false;
    }
    const QString normalized = newName.trimmed();
    if (!QRegularExpression(QStringLiteral(R"(^[A-Za-z_][A-Za-z0-9_.-]*$)"))
             .match(normalized).hasMatch()) {
        if (error) *error = QStringLiteral("patch 名称无效");
        return false;
    }
    int oldOpen = -1;
    if (!patchRange(lines_, oldName, &oldOpen, nullptr)) {
        if (error) *error = QStringLiteral("未找到 patch %1").arg(oldName);
        return false;
    }
    if (oldName != normalized && patchRange(lines_, normalized, nullptr, nullptr)) {
        if (error) *error = QStringLiteral("patch %1 已存在").arg(normalized);
        return false;
    }
    int nameLine = oldOpen;
    while (nameLine >= 0 && !stripFoamComments(lines_.at(nameLine)).contains(oldName)) --nameLine;
    if (nameLine < 0) {
        if (error) *error = QStringLiteral("patch 名称行无效");
        return false;
    }
    lines_[nameLine].replace(oldName, normalized);
    dirty_ = true;
    parse();
    return true;
}

bool ConfigDocument::removeBlockMeshPatch(const QString& name, QString* error)
{
    if (!isBlockMeshDict(path_)) {
        if (error) *error = QStringLiteral("只有 blockMeshDict 支持 patch 编辑");
        return false;
    }
    int open = -1;
    int close = -1;
    if (!patchRange(lines_, name, &open, &close)) {
        if (error) *error = QStringLiteral("未找到 patch %1").arg(name);
        return false;
    }
    int nameLine = open;
    while (nameLine >= 0 && !stripFoamComments(lines_.at(nameLine)).contains(name)) {
        --nameLine;
    }
    if (nameLine < 0) {
        if (error) *error = QStringLiteral("patch %1 的名称行无效").arg(name);
        return false;
    }
    for (int line = close; line >= nameLine; --line) lines_.removeAt(line);
    dirty_ = true;
    parse();
    return true;
}

bool ConfigDocument::removeEntry(int entryIndex, QString* error)
{
    if (entryIndex < 0 || entryIndex >= entries_.size()) {
        if (error) *error = QStringLiteral("无效的配置项索引");
        return false;
    }
    const ConfigEntry entry = entries_.at(entryIndex);
    if (entry.startLine < 0 || entry.endLine >= lines_.size()) {
        if (error) *error = QStringLiteral("配置项行范围无效");
        return false;
    }
    for (int line = entry.endLine; line >= entry.startLine; --line) lines_.removeAt(line);
    dirty_ = true;
    parse();
    return true;
}

QString ConfigDocument::fileName() const
{
    return QFileInfo(path_).fileName();
}

int ConfigDocument::squareBracketBalance(const QString& text)
{
    bool quoted = false;
    QChar quote;
    int balance = 0;
    for (int i = 0; i < text.size(); ++i) {
        const QChar ch = text.at(i);
        if ((ch == QLatin1Char('"') || ch == QLatin1Char('\''))) {
            if (!quoted) {
                quoted = true;
                quote = ch;
            } else if (quote == ch && (i == 0 || text.at(i - 1) != QLatin1Char('\\'))) {
                quoted = false;
            }
            continue;
        }
        if (quoted) continue;
        if (ch == QLatin1Char('#')) break;
        if (ch == QLatin1Char('[')) ++balance;
        if (ch == QLatin1Char(']')) --balance;
    }
    return balance;
}

QString ConfigDocument::valueWithoutComment(const QString& text)
{
    bool quoted = false;
    QChar quote;
    int squareDepth = 0;
    for (int i = 0; i < text.size(); ++i) {
        const QChar ch = text.at(i);
        if (ch == QLatin1Char('"') || ch == QLatin1Char('\'')) {
            if (!quoted) {
                quoted = true;
                quote = ch;
            } else if (quote == ch && (i == 0 || text.at(i - 1) != QLatin1Char('\\'))) {
                quoted = false;
            }
        } else if (!quoted) {
            if (ch == QLatin1Char('[')) ++squareDepth;
            if (ch == QLatin1Char(']')) --squareDepth;
            if (ch == QLatin1Char('#') && squareDepth == 0) {
                return text.left(i).trimmed();
            }
        }
    }
    return text.trimmed();
}

QString ConfigDocument::inlineComment(const QString& text)
{
    bool quoted = false;
    QChar quote;
    int squareDepth = 0;
    for (int i = 0; i < text.size(); ++i) {
        const QChar ch = text.at(i);
        if (ch == QLatin1Char('"') || ch == QLatin1Char('\'')) {
            if (!quoted) {
                quoted = true;
                quote = ch;
            } else if (quote == ch && (i == 0 || text.at(i - 1) != QLatin1Char('\\'))) {
                quoted = false;
            }
        } else if (!quoted) {
            if (ch == QLatin1Char('[')) ++squareDepth;
            if (ch == QLatin1Char(']')) --squareDepth;
            if (ch == QLatin1Char('#') && squareDepth == 0) {
                return text.mid(i).trimmed();
            }
        }
    }
    return {};
}

void ConfigDocument::parse()
{
    entries_.clear();

    if (QFileInfo(path_).fileName().compare(
            QStringLiteral("blockMeshDict"), Qt::CaseInsensitive) == 0) {
        parseBlockMeshDict();
        return;
    }

    if (lines_.join(QLatin1Char('\n')).contains(
            QRegularExpression(QStringLiteral(R"(\bFoamFile\b)")))) {
        parseFoamDictionary();
        return;
    }

    QHash<QString, int> sectionCounts;
    QString section;
    QString sectionLabel;
    bool caseMode = false;

    const QRegularExpression sectionPattern(
        QStringLiteral(R"(^\s*\[([^\]]+)\]\s*(?:#.*)?$)"));
    const QRegularExpression assignmentPattern(
        QStringLiteral(R"(^\s*([A-Za-z_][A-Za-z0-9_.-]*)\s*=\s*(.*)$)"));

    for (int i = 0; i < lines_.size(); ++i) {
        const QString line = lines_.at(i);
        const QString trimmed = line.trimmed();
        if (trimmed.isEmpty()) continue;

        if (trimmed.startsWith(QLatin1Char('#'))) {
            caseMode = true;
            section = trimmed.mid(1).trimmed();
            const int occurrence = ++sectionCounts[section];
            sectionLabel = occurrence > 1
                ? QStringLiteral("%1 %2").arg(section).arg(occurrence)
                : section;

            int valueLine = i + 1;
            while (valueLine < lines_.size()
                   && lines_.at(valueLine).trimmed().isEmpty()) {
                ++valueLine;
            }
            if (valueLine < lines_.size()
                && !lines_.at(valueLine).trimmed().startsWith(QLatin1Char('#'))) {
                ConfigEntry entry;
                entry.section = section;
                entry.sectionLabel = sectionLabel;
                entry.key = QStringLiteral("value");
                entry.value = lines_.at(valueLine).trimmed();
                entry.startLine = valueLine;
                entry.endLine = valueLine;
                entry.caseSection = true;
                entries_.push_back(entry);
                i = valueLine;
            }
            continue;
        }

        const QRegularExpressionMatch sectionMatch = sectionPattern.match(line);
        if (sectionMatch.hasMatch()) {
            caseMode = false;
            section = sectionMatch.captured(1).trimmed();
            const int occurrence = ++sectionCounts[section];
            sectionLabel = occurrence > 1
                ? QStringLiteral("%1 %2").arg(section).arg(occurrence)
                : section;
            continue;
        }

        if (caseMode || section.isEmpty()) continue;

        const QRegularExpressionMatch assignmentMatch =
            assignmentPattern.match(line);
        if (!assignmentMatch.hasMatch()) continue;

        ConfigEntry entry;
        entry.section = section;
        entry.sectionLabel = sectionLabel;
        entry.key = assignmentMatch.captured(1);
        entry.startLine = i;
        entry.endLine = i;

        QString value = valueWithoutComment(assignmentMatch.captured(2));
        int balance = squareBracketBalance(value);
        while (balance > 0 && entry.endLine + 1 < lines_.size()) {
            ++entry.endLine;
            const QString continuation = lines_.at(entry.endLine);
            value += QLatin1Char('\n') + continuation;
            balance += squareBracketBalance(continuation);
        }
        entry.value = value.trimmed();
        entries_.push_back(entry);
        i = entry.endLine;
    }
}

void ConfigDocument::parseFoamDictionary()
{
    QStringList blocks;
    QString pendingBlock;
    const QRegularExpression assignmentPattern(
        QStringLiteral(R"(^\s*([A-Za-z_][A-Za-z0-9_.-]*)\s+(.+?)\s*;\s*$)"));
    const QRegularExpression blockNamePattern(
        QStringLiteral(R"(^\s*([A-Za-z_][A-Za-z0-9_.-]*)\s*$)"));

    for (int i = 0; i < lines_.size(); ++i) {
        const QString clean = stripFoamComments(lines_.at(i)).trimmed();
        if (clean.isEmpty()) continue;

        if (clean.startsWith(QLatin1Char('}'))) {
            if (!blocks.isEmpty()) blocks.removeLast();
            pendingBlock.clear();
            continue;
        }

        const int brace = clean.indexOf(QLatin1Char('{'));
        if (brace >= 0) {
            QString blockName = clean.left(brace).trimmed();
            if (blockName.isEmpty()) blockName = pendingBlock;
            if (!blockName.isEmpty()) blocks.push_back(blockName);
            pendingBlock.clear();
            continue;
        }

        const QRegularExpressionMatch assignment = assignmentPattern.match(clean);
        if (assignment.hasMatch()) {
            // FoamFile 只是 OpenFOAM 文件头元数据，不属于用户可编辑的业务配置。
            if (!blocks.isEmpty()
                && blocks.constFirst().compare(QStringLiteral("FoamFile"),
                                               Qt::CaseInsensitive) == 0) {
                pendingBlock.clear();
                continue;
            }
            ConfigEntry entry;
            entry.section = blocks.isEmpty()
                ? QString() : blocks.join(QLatin1Char('/'));
            entry.sectionLabel = entry.section;
            entry.key = assignment.captured(1);
            entry.value = assignment.captured(2).trimmed();
            entry.startLine = i;
            entry.endLine = i;
            entry.foamAssignment = true;
            entries_.push_back(entry);
            pendingBlock.clear();
            continue;
        }

        const QRegularExpressionMatch blockName = blockNamePattern.match(clean);
        pendingBlock = blockName.hasMatch() ? blockName.captured(1) : QString();
    }
}

void ConfigDocument::parseBlockMeshDict()
{
    const QRegularExpression vertexLine(
        QStringLiteral(R"(^\s*\(\s*[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s*\)\s*$)"));
    const QRegularExpression blockLine(
        QStringLiteral(R"(^\s*(hex|wedge|prism)\s*\([^\)]*\).*$)"),
        QRegularExpression::CaseInsensitiveOption);
    const QRegularExpression edgeLine(
        QStringLiteral(R"(^\s*(arc|spline|polyLine|polySpline|line)\b.*$)"),
        QRegularExpression::CaseInsensitiveOption);
    const QRegularExpression faceLine(
        QStringLiteral(R"(^\s*\(\s*\d+(?:\s+\d+)+\s*\)\s*$)"));
    const QRegularExpression typeLine(
        QStringLiteral(R"(^\s*type\s+([^;]+);\s*$)"),
        QRegularExpression::CaseInsensitiveOption);

    int vertexIndex = 0;
    int blockIndex = 0;
    int edgeIndex = 0;

    for (int i = 0; i < lines_.size(); ++i) {
        const QString trimmed = stripFoamComments(lines_.at(i)).trimmed();
        if (trimmed.compare(QStringLiteral("vertices"),
                            Qt::CaseInsensitive) == 0
            || trimmed.startsWith(QStringLiteral("vertices("),
                                  Qt::CaseInsensitive)) {
            const int end = findMatching(lines_, i,
                                         QLatin1Char('('), QLatin1Char(')'));
            for (int line = i + 1; line >= 0 && line < end; ++line) {
                const QString value = stripFoamComments(lines_.at(line)).trimmed();
                if (!vertexLine.match(value).hasMatch()) continue;
                addLineEntry(&entries_, QStringLiteral("vertices"),
                             QStringLiteral("v%1").arg(vertexIndex++),
                             value, line);
                entries_.last().coordinateTuple = true;
            }
            i = std::max(i, end);
            continue;
        }

        if (trimmed.compare(QStringLiteral("blocks"),
                            Qt::CaseInsensitive) == 0
            || trimmed.startsWith(QStringLiteral("blocks("),
                                  Qt::CaseInsensitive)) {
            const int end = findMatching(lines_, i,
                                         QLatin1Char('('), QLatin1Char(')'));
            for (int line = i + 1; line >= 0 && line < end; ++line) {
                const QString value = stripFoamComments(lines_.at(line)).trimmed();
                if (!blockLine.match(value).hasMatch()) continue;
                addLineEntry(&entries_, QStringLiteral("blocks"),
                             QStringLiteral("block%1").arg(blockIndex++),
                             value, line);
            }
            i = std::max(i, end);
            continue;
        }

        if (trimmed.compare(QStringLiteral("edges"),
                            Qt::CaseInsensitive) == 0
            || trimmed.startsWith(QStringLiteral("edges("),
                                  Qt::CaseInsensitive)) {
            const int end = findMatching(lines_, i,
                                         QLatin1Char('('), QLatin1Char(')'));
            for (int line = i + 1; line >= 0 && line < end; ++line) {
                const QString value = stripFoamComments(lines_.at(line)).trimmed();
                if (!edgeLine.match(value).hasMatch()) continue;
                addLineEntry(&entries_, QStringLiteral("edges"),
                             QStringLiteral("edge%1").arg(edgeIndex++),
                             value, line);
            }
            i = std::max(i, end);
            continue;
        }

        if (trimmed.compare(QStringLiteral("boundary"),
                            Qt::CaseInsensitive) != 0
            && !trimmed.startsWith(QStringLiteral("boundary("),
                                   Qt::CaseInsensitive)) {
            continue;
        }

        const int boundaryEnd = findMatching(lines_, i,
                                             QLatin1Char('('),
                                             QLatin1Char(')'));
        int line = i + 1;
        while (line > 0 && line < boundaryEnd) {
            QString patchName = stripFoamComments(lines_.at(line)).trimmed();
            if (patchName.isEmpty() || patchName == QStringLiteral("(")
                || patchName == QStringLiteral(")")) {
                ++line;
                continue;
            }
            if (patchName.contains(QLatin1Char('{'))) {
                patchName = patchName.left(patchName.indexOf(QLatin1Char('{')))
                                .trimmed();
            }

            int openLine = line;
            if (!stripFoamComments(lines_.at(openLine)).contains(QLatin1Char('{'))) {
                ++openLine;
            }
            if (openLine >= boundaryEnd
                || !stripFoamComments(lines_.at(openLine)).contains(QLatin1Char('{'))) {
                ++line;
                continue;
            }
            const int patchEnd = findMatching(lines_, openLine,
                                              QLatin1Char('{'),
                                              QLatin1Char('}'));
            if (patchEnd < 0) break;

            int faceIndex = 0;
            for (int patchLine = openLine + 1;
                 patchLine < patchEnd; ++patchLine) {
                const QString value =
                    stripFoamComments(lines_.at(patchLine)).trimmed();
                const QRegularExpressionMatch type = typeLine.match(value);
                if (type.hasMatch()) {
                    addLineEntry(&entries_,
                                 QStringLiteral("boundary/%1").arg(patchName),
                                 QStringLiteral("type"),
                                 value, patchLine);
                } else if (faceLine.match(value).hasMatch()) {
                    addLineEntry(&entries_,
                                 QStringLiteral("boundary/%1").arg(patchName),
                                 QStringLiteral("face%1").arg(faceIndex++),
                                 value, patchLine);
                }
            }
            line = patchEnd + 1;
        }
        i = std::max(i, boundaryEnd);
    }
}

} // namespace SF::GUI
