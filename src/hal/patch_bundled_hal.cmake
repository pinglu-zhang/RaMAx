if(NOT DEFINED RAMAX_HAL_SOURCE_DIR OR
   NOT DEFINED RAMAX_HAL_PATCH_OUTPUT_DIR)
    message(FATAL_ERROR "HAL patch paths were not configured")
endif()

function(_ramax_read_pinned_hal_source source expected_sha output_variable)
    file(SHA256 "${source}" actual_sha)
    if(NOT actual_sha STREQUAL expected_sha)
        message(FATAL_ERROR
            "Bundled HAL source drifted at ${source}\n"
            "expected SHA256 ${expected_sha}\n"
            "actual   SHA256 ${actual_sha}\n"
            "Review and update the RaMAx safety patch explicitly.")
    endif()
    file(READ "${source}" contents)
    set(${output_variable} "${contents}" PARENT_SCOPE)
endfunction()

function(_ramax_copy_pinned_hal_source relative_path expected_sha)
    set(source "${RAMAX_HAL_SOURCE_DIR}/${relative_path}")
    _ramax_read_pinned_hal_source("${source}" "${expected_sha}" contents)
    get_filename_component(output_name "${relative_path}" NAME)
    file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/${output_name}" "${contents}")
endfunction()

macro(_ramax_replace_hal_block contents_variable before_variable after_variable label)
    string(FIND "${${contents_variable}}" "${${before_variable}}" match_offset)
    if(match_offset EQUAL -1)
        message(FATAL_ERROR "Pinned HAL patch block not found: ${label}")
    endif()
    string(LENGTH "${${before_variable}}" match_length)
    math(EXPR remaining_offset "${match_offset} + ${match_length}")
    string(SUBSTRING "${${contents_variable}}" ${remaining_offset} -1 remaining)
    string(FIND "${remaining}" "${${before_variable}}" duplicate_offset)
    if(NOT duplicate_offset EQUAL -1)
        message(FATAL_ERROR "Pinned HAL patch block is ambiguous: ${label}")
    endif()
    string(REPLACE "${${before_variable}}" "${${after_variable}}"
        ${contents_variable} "${${contents_variable}}")
endmacro()

set(alignment_source
    "${RAMAX_HAL_SOURCE_DIR}/api/hdf5_impl/hdf5Alignment.cpp")
_ramax_read_pinned_hal_source(
    "${alignment_source}"
    "43b15f415b45b38ffb8291a0be6af46d2c60bce9690325257d945ad82709521e"
    alignment_contents)

set(before_includes [=[#include <deque>
#include <fstream>
#include <iostream>]=])
set(after_includes [=[#include <deque>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>]=])
_ramax_replace_hal_block(alignment_contents before_includes after_includes
    "alignment exception support")

set(before_constructors [=[Hdf5Alignment::Hdf5Alignment(const string &alignmentPath, unsigned mode, const H5::FileCreatPropList &fileCreateProps,
                             const H5::FileAccPropList &fileAccessProps, const H5::DSetCreatPropList &datasetCreateProps,
                             bool inMemory)
    : _alignmentPath(alignmentPath), _mode(halDefaultAccessMode(mode)), _file(NULL), _flags(hdf5DefaultFlags(_mode)),
      _inMemory(inMemory), _metaData(NULL), _tree(NULL), _dirty(false) {
    _cprops.copy(fileCreateProps);
    _aprops.copy(fileAccessProps);
    _dcprops.copy(datasetCreateProps);
    if (_inMemory) {
        setInMemory();
    }
    if (_mode & CREATE_ACCESS) {
        create();
    } else {
        open();
    }
}

Hdf5Alignment::Hdf5Alignment(const std::string &alignmentPath, unsigned mode, const CLParser *parser)
    : _alignmentPath(alignmentPath), _mode(halDefaultAccessMode(mode)), _file(NULL), _flags(hdf5DefaultFlags(_mode)),
      _inMemory(false), _metaData(NULL), _tree(NULL), _dirty(false) {
    initializeFromOptions(parser);
    if (_inMemory) {
        setInMemory();
    }
    if (_mode & CREATE_ACCESS) {
        create();
    } else {
        open();
    }
}]=])
set(after_constructors [=[Hdf5Alignment::Hdf5Alignment(const string &alignmentPath, unsigned mode, const H5::FileCreatPropList &fileCreateProps,
                             const H5::FileAccPropList &fileAccessProps, const H5::DSetCreatPropList &datasetCreateProps,
                             bool inMemory)
    : _alignmentPath(alignmentPath), _mode(halDefaultAccessMode(mode)), _file(NULL), _flags(hdf5DefaultFlags(_mode)),
      _inMemory(inMemory), _metaData(NULL), _tree(NULL), _dirty(false) {
    try {
        _cprops.copy(fileCreateProps);
        _aprops.copy(fileAccessProps);
        _dcprops.copy(datasetCreateProps);
        if (_inMemory) {
            setInMemory();
        }
        if (_mode & CREATE_ACCESS) {
            create();
        } else {
            open();
        }
    } catch (...) {
        exception_ptr originalFailure = current_exception();
        try {
            close();
        } catch (...) {
        }
        rethrow_exception(originalFailure);
    }
}

Hdf5Alignment::Hdf5Alignment(const std::string &alignmentPath, unsigned mode, const CLParser *parser)
    : _alignmentPath(alignmentPath), _mode(halDefaultAccessMode(mode)), _file(NULL), _flags(hdf5DefaultFlags(_mode)),
      _inMemory(false), _metaData(NULL), _tree(NULL), _dirty(false) {
    try {
        initializeFromOptions(parser);
        if (_inMemory) {
            setInMemory();
        }
        if (_mode & CREATE_ACCESS) {
            create();
        } else {
            open();
        }
    } catch (...) {
        exception_ptr originalFailure = current_exception();
        try {
            close();
        } catch (...) {
        }
        rethrow_exception(originalFailure);
    }
}]=])
_ramax_replace_hal_block(alignment_contents before_constructors after_constructors
    "alignment constructor cleanup")

set(before_destructor [=[Hdf5Alignment::~Hdf5Alignment() {
    close();
}]=])
set(after_destructor [=[Hdf5Alignment::~Hdf5Alignment() {
    try {
        close();
    } catch (...) {
    }
}]=])
_ramax_replace_hal_block(alignment_contents before_destructor after_destructor
    "alignment non-throwing destruction")

set(before_close [=[void Hdf5Alignment::close() {
    if (_file != NULL) {
        if (not isReadOnly()) {
            writeTree();
        }
        if (_tree != NULL) {
            stTree_destruct(_tree);
            _tree = NULL;
        }
        // todo: make sure there's no memory leak with metadata
        // smart pointer should prevent
        if (_metaData != NULL) {
            if (not isReadOnly()) {
                _metaData->write();
            }
            delete _metaData;
            _metaData = NULL;
        }
        writeVersion();
        map<string, Hdf5Genome *>::iterator mapIt;
        for (mapIt = _openGenomes.begin(); mapIt != _openGenomes.end(); ++mapIt) {
            Hdf5Genome *genome = mapIt->second;
            if (not isReadOnly()) {
                genome->write();
            }
            delete genome;
        }
        _openGenomes.clear();
        if (not isReadOnly()) {
            _file->flush(H5F_SCOPE_LOCAL);
        }
        _file->close();
        delete _file;
        _file = NULL;
    } else {
        assert(_tree == NULL);
        assert(_openGenomes.empty() == true);
    }
}]=])
set(after_close [=[void Hdf5Alignment::close() {
    exception_ptr firstFailure;
    const bool readOnly = isReadOnly();
    const auto retainFailure = [&firstFailure](const auto &operation) {
        try {
            operation();
        } catch (...) {
            if (not firstFailure) {
                firstFailure = current_exception();
            }
        }
    };

    if (_file != NULL && not readOnly) {
        retainFailure([this]() { writeTree(); });
    }

    if (_tree != NULL) {
        stTree_destruct(_tree);
        _tree = NULL;
    }
    _nodeMap.clear();

    HDF5MetaData *metaData = _metaData;
    _metaData = NULL;
    if (metaData != NULL) {
        if (not readOnly) {
            retainFailure([metaData]() { metaData->write(); });
        }
        delete metaData;
    }

    if (_file != NULL) {
        retainFailure([this]() { writeVersion(); });
    }

    while (not _openGenomes.empty()) {
        map<string, Hdf5Genome *>::iterator mapIt = _openGenomes.begin();
        Hdf5Genome *genome = mapIt->second;
        _openGenomes.erase(mapIt);
        if (not readOnly) {
            retainFailure([genome]() { genome->write(); });
        }
        delete genome;
    }

    H5File *file = _file;
    _file = NULL;
    if (file != NULL) {
        if (not readOnly) {
            retainFailure([file]() { file->flush(H5F_SCOPE_LOCAL); });
        }
        retainFailure([file]() { file->close(); });
        delete file;
    }
    _dirty = false;

    if (firstFailure) {
        rethrow_exception(firstFailure);
    }
}]=])
_ramax_replace_hal_block(alignment_contents before_close after_close
    "alignment complete close cleanup")

set(before_close_genome [=[void Hdf5Alignment::closeGenome(const Genome *genome) const {
    string name = genome->getName();
    map<string, Hdf5Genome *>::iterator mapIt = _openGenomes.find(name);
    if (mapIt == _openGenomes.end()) {
        throw hal_exception("Attempt to close non-open genome.  "
                            "Should not even be possible");
    }
    mapIt->second->write();
    delete mapIt->second;
    _openGenomes.erase(mapIt);

    // reset the parent/child genoem cachces (which store genome pointers to
    // the genome we're closing
    if (name != getRootName()) {
        mapIt = _openGenomes.find(getParentName(name));
        if (mapIt != _openGenomes.end()) {
            mapIt->second->resetBranchCaches();
        }
    }
    vector<string> childNames = getChildNames(name);
    for (size_t i = 0; i < childNames.size(); ++i) {
        mapIt = _openGenomes.find(childNames[i]);
        if (mapIt != _openGenomes.end()) {
            mapIt->second->resetBranchCaches();
        }
    }
}]=])
set(after_close_genome [=[void Hdf5Alignment::closeGenome(const Genome *genome) const {
    string name = genome->getName();
    map<string, Hdf5Genome *>::iterator mapIt = _openGenomes.find(name);
    if (mapIt == _openGenomes.end()) {
        throw hal_exception("Attempt to close non-open genome.  "
                            "Should not even be possible");
    }

    // Resolve allocating tree lookups before consuming the caller's handle.
    // A failed lookup or write leaves the genome open and safe to close again.
    const bool hasParent = name != getRootName();
    const string parentName = hasParent ? getParentName(name) : string();
    const vector<string> childNames = getChildNames(name);
    mapIt->second->write();
    delete mapIt->second;
    _openGenomes.erase(mapIt);

    if (hasParent) {
        mapIt = _openGenomes.find(parentName);
        if (mapIt != _openGenomes.end()) {
            mapIt->second->resetBranchCaches();
        }
    }
    for (const string &childName : childNames) {
        mapIt = _openGenomes.find(childName);
        if (mapIt != _openGenomes.end()) {
            mapIt->second->resetBranchCaches();
        }
    }
}]=])
_ramax_replace_hal_block(alignment_contents before_close_genome after_close_genome
    "genome close preserves ownership on failure")

set(before_tree_write [=[void Hdf5Alignment::writeTree() {
    if (_dirty == false)
        return;

    char *treeString = NULL;
    if (_tree != NULL) {
        treeString = stTree_getNewickTreeString(_tree);
    } else {
        treeString = (char *)malloc(sizeof(char));
        treeString[0] = '\0';
    }
    assert(_file != NULL);
    HDF5MetaData treeMeta(_file, TreeGroupName);
    treeMeta.set(TreeGroupName, treeString);
    free(treeString);
}]=])
set(after_tree_write [=[void Hdf5Alignment::writeTree() {
    if (_dirty == false)
        return;

    unique_ptr<char, decltype(&std::free)> treeString(NULL, &std::free);
    if (_tree != NULL) {
        treeString.reset(stTree_getNewickTreeString(_tree));
    }
    assert(_file != NULL);
    HDF5MetaData treeMeta(_file, TreeGroupName);
    treeMeta.set(TreeGroupName, treeString ? treeString.get() : "");
    treeMeta.write();
}]=])
_ramax_replace_hal_block(alignment_contents before_tree_write after_tree_write
    "exception-safe explicit tree metadata write")

set(before_version_write [=[    HDF5MetaData versionMeta(_file, VersionGroupName);
    versionMeta.set(VersionGroupName, HAL_VERSION);]=])
set(after_version_write [=[    HDF5MetaData versionMeta(_file, VersionGroupName);
    versionMeta.set(VersionGroupName, HAL_VERSION);
    versionMeta.write();]=])
_ramax_replace_hal_block(alignment_contents before_version_write after_version_write
    "explicit version metadata write")

set(genome_source
    "${RAMAX_HAL_SOURCE_DIR}/api/hdf5_impl/hdf5Genome.cpp")
_ramax_read_pinned_hal_source(
    "${genome_source}"
    "619345c30a527c21828588cb4a77d6e4a33d583fe3550ff8b80d7d2add1db9d9"
    genome_contents)
set(before_genome_includes [=[#include "hdf5Genome.h"
#include "H5Cpp.h"]=])
set(after_genome_includes [=[#include "hdf5Genome.h"
#include "ramaxHdf5BulkDna.h"
#include "H5Cpp.h"]=])
_ramax_replace_hal_block(genome_contents before_genome_includes
    after_genome_includes "RaMAx fused segment writer declaration")
set(before_genome_standard_includes [=[#include <algorithm>
#include <cassert>
#include <iostream>]=])
set(after_genome_standard_includes [=[#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>]=])
_ramax_replace_hal_block(genome_contents before_genome_standard_includes
    after_genome_standard_includes "RaMAx fused field copy support")
set(before_genome_constructor [=[Hdf5Genome::Hdf5Genome(const string &name, Hdf5Alignment *alignment, PortableH5Location *h5Parent,
                       const DSetCreatPropList &dcProps, bool inMemory)
    : Genome(alignment, name), _alignment(alignment), _h5Parent(h5Parent), _name(name), _numChildrenInBottomArray(0),
      _totalSequenceLength(0), _numChunksInArrayBuffer(inMemory ? 0 : 1) {
    _dcprops.copy(dcProps);
    assert(!name.empty());
    assert(alignment != NULL && h5Parent != NULL);

    try {
        HDF5DisableExceptionPrinting prDisable;
        _group = h5Parent->openGroup(name);
    } catch (Exception &e) {
        _group = h5Parent->createGroup(name);
    }
    read();
    _metaData = new HDF5MetaData(&_group, metaGroupName);
    _rup = new HDF5MetaData(&_group, rupGroupName);

    _totalSequenceLength = _dnaArray.getSize() * 2;
    if (_totalSequenceLength > 0 && _rup->get(rupGroupName) == "1") {
        _totalSequenceLength -= 1;
    } else if (_totalSequenceLength == 0 && _sequenceIdxArray.getSize() > 0) {
        Hdf5Sequence lastSeq(this, &_sequenceIdxArray, &_sequenceNameArray, _sequenceNameArray.getSize() - 1);
        _totalSequenceLength = lastSeq.getEndPosition() + 1;
    }
}]=])
set(after_genome_constructor [=[Hdf5Genome::Hdf5Genome(const string &name, Hdf5Alignment *alignment, PortableH5Location *h5Parent,
                       const DSetCreatPropList &dcProps, bool inMemory)
    : Genome(alignment, name), _alignment(alignment), _h5Parent(h5Parent), _name(name), _metaData(NULL), _rup(NULL),
      _numChildrenInBottomArray(0), _totalSequenceLength(0), _numChunksInArrayBuffer(inMemory ? 0 : 1) {
    try {
        _dcprops.copy(dcProps);
        assert(!name.empty());
        assert(alignment != NULL && h5Parent != NULL);

        try {
            HDF5DisableExceptionPrinting prDisable;
            _group = h5Parent->openGroup(name);
        } catch (Exception &e) {
            _group = h5Parent->createGroup(name);
        }
        read();
        _metaData = new HDF5MetaData(&_group, metaGroupName);
        _rup = new HDF5MetaData(&_group, rupGroupName);

        _totalSequenceLength = _dnaArray.getSize() * 2;
        if (_totalSequenceLength > 0 && _rup->get(rupGroupName) == "1") {
            _totalSequenceLength -= 1;
        } else if (_totalSequenceLength == 0 && _sequenceIdxArray.getSize() > 0) {
            Hdf5Sequence lastSeq(this, &_sequenceIdxArray, &_sequenceNameArray, _sequenceNameArray.getSize() - 1);
            _totalSequenceLength = lastSeq.getEndPosition() + 1;
        }
    } catch (...) {
        delete _rup;
        _rup = NULL;
        delete _metaData;
        _metaData = NULL;
        throw;
    }
}]=])
_ramax_replace_hal_block(genome_contents before_genome_constructor
    after_genome_constructor "genome constructor cleanup")

set(genome_header
    "${RAMAX_HAL_SOURCE_DIR}/api/hdf5_impl/hdf5Genome.h")
_ramax_read_pinned_hal_source(
    "${genome_header}"
    "fb6a93537db39cccdfee9f80d27c37dd9297c6c3619c4dfbffa331cd4f8ec6d5"
    genome_header_contents)
set(before_genome_forward_declarations [=[    class Hdf5Sequence;]=])
set(after_genome_forward_declarations [=[    class Hdf5Sequence;
    class RamaxHdf5SegmentWriter;
    void ramaxWriteBulkDna(Genome &genome, hal_index_t start,
                           const char *dna, hal_size_t length);]=])
_ramax_replace_hal_block(genome_header_contents
    before_genome_forward_declarations after_genome_forward_declarations
    "RaMAx bulk DNA namespace declaration")
set(before_genome_friends [=[        friend class Hdf5TopSegment;
        friend class Hdf5BottomSegment;
        friend class Hdf5SequenceIterator;
        friend class Hdf5Sequence;]=])
set(after_genome_friends [=[        friend class Hdf5TopSegment;
        friend class Hdf5BottomSegment;
        friend class Hdf5SequenceIterator;
        friend class Hdf5Sequence;
        friend class RamaxHdf5SegmentWriter;
        friend void ramaxWriteBulkDna(Genome &genome, hal_index_t start,
                                      const char *dna, hal_size_t length);]=])
_ramax_replace_hal_block(genome_header_contents before_genome_friends
    after_genome_friends "RaMAx bounded bulk DNA access")

set(bottom_segment_header
    "${RAMAX_HAL_SOURCE_DIR}/api/hdf5_impl/hdf5BottomSegment.h")
_ramax_read_pinned_hal_source(
    "${bottom_segment_header}"
    "b6574f7116c1798e01730ed4494e16cd40339ce778063b472f156043753aa08e"
    bottom_segment_header_contents)
set(before_bottom_segment_private [=[      private:
        Hdf5Genome *getHdf5Genome() const {]=])
set(after_bottom_segment_private [=[      private:
        friend class RamaxHdf5SegmentWriter;

        Hdf5Genome *getHdf5Genome() const {]=])
_ramax_replace_hal_block(bottom_segment_header_contents
    before_bottom_segment_private after_bottom_segment_private
    "RaMAx fused bottom-row access")

set(top_segment_header
    "${RAMAX_HAL_SOURCE_DIR}/api/hdf5_impl/hdf5TopSegment.h")
_ramax_read_pinned_hal_source(
    "${top_segment_header}"
    "29b451d10f75cc2c4368de75bd0ec03b5e6878796f8fe921744b1161e390e170"
    top_segment_header_contents)
set(before_top_segment_private [=[      private:
        Hdf5Genome *getHdf5Genome() const {]=])
set(after_top_segment_private [=[      private:
        friend class RamaxHdf5SegmentWriter;

        Hdf5Genome *getHdf5Genome() const {]=])
_ramax_replace_hal_block(top_segment_header_contents
    before_top_segment_private after_top_segment_private
    "RaMAx fused top-row access")

set(dna_driver_header
    "${RAMAX_HAL_SOURCE_DIR}/api/hdf5_impl/hdf5DnaDriver.h")
_ramax_read_pinned_hal_source(
    "${dna_driver_header}"
    "43b846ea3d5f4bad75891d302719678ae4ea0868559a7a5b3b66dac5ba2bd36d"
    dna_driver_header_contents)
set(before_dna_driver_flush [=[        void flush();]=])
set(after_dna_driver_flush [=[        void flush();

        // RaMAx-owned bulk path. It writes directly into the active packed
        // page while retaining DnaAccess's paging and dirtiness protocol.
        void ramaxWriteString(hal_index_t start, const char *dna,
                              hal_size_t length);]=])
_ramax_replace_hal_block(dna_driver_header_contents before_dna_driver_flush
    after_dna_driver_flush "bounded packed DNA writer declaration")

set(dna_driver_source
    "${RAMAX_HAL_SOURCE_DIR}/api/hdf5_impl/hdf5DnaDriver.cpp")
_ramax_read_pinned_hal_source(
    "${dna_driver_source}"
    "b4fd8729376a78b08806766448e8b0193b34688f86d68f69073f712289eff052"
    dna_driver_contents)
string(PREPEND dna_driver_contents "#include <array>\n")
string(APPEND dna_driver_contents [=[

namespace {
const std::array<uint8_t, 256> ramaxPackedDnaCodes = [] {
    std::array<uint8_t, 256> codes;
    codes.fill(0xFFU);
    for (const unsigned char base : {'A', 'C', 'G', 'T', 'N',
                                    'a', 'c', 'g', 't', 'n'}) {
        codes[base] = hal::dnaPackMap[base];
    }
    for (const unsigned char base : {'K', 'M', 'R', 'Y', 'U', 'S', 'W',
                                    'B', 'D', 'H', 'V'}) {
        codes[base] = hal::dnaPackMap[static_cast<uint8_t>('N')];
        codes[base - 'A' + 'a'] = hal::dnaPackMap[static_cast<uint8_t>('n')];
    }
    return codes;
}();

uint8_t ramaxPackedDnaCode(char base) {
    const uint8_t code = ramaxPackedDnaCodes[static_cast<uint8_t>(base)];
    if (code == 0xFFU) {
        throw ::hal_exception(std::string("Trying to set invalid character: ") + base);
    }
    return code;
}
}

void HDF5DnaAccess::ramaxWriteString(hal_index_t start, const char *dna,
                                    hal_size_t length) {
    hal_size_t inputOffset = 0;
    hal_index_t index = start;
    // HAL's format is byte-order independent: even genome positions occupy
    // the high nibble and odd positions the low nibble. Pages therefore start
    // on an even position. Partial first/last bytes retain the neighboring
    // nibble, including the unused odd-genome-length sentinel nibble.
    while (inputOffset < length) {
        const hal_index_t relativeStart = access(index);
        const hal_size_t pageBases =
            static_cast<hal_size_t>(_endIndex - index);
        const hal_size_t remaining = length - inputOffset;
        hal_size_t count = remaining < pageBases ? remaining : pageBases;
        unsigned char *output = reinterpret_cast<unsigned char *>(
            _buffer + relativeStart / 2);

        if ((relativeStart & 1) != 0 && count != 0) {
            const uint8_t code = ramaxPackedDnaCode(dna[inputOffset]);
            _dirty = true;
            *output = static_cast<unsigned char>((*output & 0xF0U) | code);
            ++output;
            ++inputOffset;
            ++index;
            --count;
        }
        while (count >= 2) {
            const uint8_t high = ramaxPackedDnaCode(dna[inputOffset]);
            const uint8_t low = ramaxPackedDnaCode(dna[inputOffset + 1]);
            _dirty = true;
            *output++ = static_cast<unsigned char>((high << 4) | low);
            inputOffset += 2;
            index += 2;
            count -= 2;
        }
        if (count != 0) {
            const uint8_t code = ramaxPackedDnaCode(dna[inputOffset]);
            _dirty = true;
            *output = static_cast<unsigned char>((*output & 0x0FU) |
                                                 (code << 4));
            ++inputOffset;
            ++index;
        }
    }
}
]=])

string(APPEND genome_contents [=[

void hal::ramaxWriteBulkDna(Genome &genome, hal_index_t start,
                            const char *dna, hal_size_t length) {
    Hdf5Genome *hdf5Genome = dynamic_cast<Hdf5Genome *>(&genome);
    if (hdf5Genome == NULL) {
        throw hal_exception("RaMAx bulk DNA writer requires an HDF5 genome");
    }
    if (start < 0 ||
        static_cast<hal_size_t>(start) > hdf5Genome->_totalSequenceLength ||
        length > hdf5Genome->_totalSequenceLength -
                     static_cast<hal_size_t>(start)) {
        throw hal_exception("RaMAx bulk DNA write is out of genome bounds");
    }
    if (length == 0) {
        return;
    }
    if (dna == NULL) {
        throw hal_exception("RaMAx bulk DNA write has a null input buffer");
    }

    HDF5DnaAccess *access =
        dynamic_cast<HDF5DnaAccess *>(hdf5Genome->_dnaAccess.get());
    if (access == NULL) {
        throw hal_exception("RaMAx bulk DNA writer found no HDF5 DNA array");
    }
    try {
        access->ramaxWriteString(start, dna, length);
        access->flush();
    } catch (...) {
        try {
            access->flush();
        } catch (...) {
        }
        throw;
    }
}

namespace {
template <typename T>
void ramaxStoreSegmentField(char *row, size_t offset,
                            const T &value) noexcept {
    std::memcpy(row + offset, &value, sizeof(value));
}

char *ramaxReacquireSegmentRow(Hdf5ExternalArray &array, hsize_t index,
                               char *row) {
    if (index < array.getBufStart() || index > array.getBufEnd()) {
        return array.getUpdate(index);
    }
    return row;
}
}

hal::RamaxHdf5SegmentWriter::RamaxHdf5SegmentWriter(Genome &genome)
    : _genome(dynamic_cast<Hdf5Genome *>(&genome)), _numChildren(0) {
    if (_genome == NULL) {
        throw ::hal_exception(
            "RaMAx fused segment writer requires an HDF5 genome");
    }
    // The on-disk child-count cache is populated by read(), not setDimensions().
    // A freshly created genome needs the live tree count, as the scalar setters do.
    _numChildren = _genome->getNumChildren();
}

void hal::RamaxHdf5SegmentWriter::initializeTop(
    hal_index_t arrayIndex, hal_index_t start, hal_size_t length,
    bool parentReversed) {
    if (start >= static_cast<hal_index_t>(_genome->_totalSequenceLength) ||
        start + length > _genome->_totalSequenceLength) {
        throw ::hal_exception(
            "Trying to set top segment coordinate out of range");
    }

    Hdf5ExternalArray &array = _genome->_topArray;
    const hsize_t rowIndex = static_cast<hsize_t>(arrayIndex);
    char *row = array.getUpdate(rowIndex);
    ramaxStoreSegmentField(
        row, Hdf5TopSegment::genomeIndexOffset, start);
    const auto end = start + length;
    ramaxStoreSegmentField(
        array.getUpdate(rowIndex + 1),
        Hdf5TopSegment::genomeIndexOffset, end);
    row = ramaxReacquireSegmentRow(array, rowIndex, row);

    const hal_index_t nullIndex = NULL_INDEX;
    ramaxStoreSegmentField(
        row, Hdf5TopSegment::parentReversedOffset, parentReversed);
    ramaxStoreSegmentField(
        row, Hdf5TopSegment::parentIndexOffset, nullIndex);
    ramaxStoreSegmentField(
        row, Hdf5TopSegment::bottomIndexOffset, nullIndex);
    ramaxStoreSegmentField(
        row, Hdf5TopSegment::parIndexOffset, nullIndex);
}

void hal::RamaxHdf5SegmentWriter::initializeBottom(
    hal_index_t arrayIndex, hal_index_t start, hal_size_t length) {
    if (start >= static_cast<hal_index_t>(_genome->_totalSequenceLength) ||
        start + length > _genome->_totalSequenceLength) {
        throw ::hal_exception(
            "Trying to set bottom segment coordinate out of range");
    }

    Hdf5ExternalArray &array = _genome->_bottomArray;
    const hsize_t rowIndex = static_cast<hsize_t>(arrayIndex);
    char *row = array.getUpdate(rowIndex);
    ramaxStoreSegmentField(
        row, Hdf5BottomSegment::genomeIndexOffset, start);
    const auto end = start + length;
    ramaxStoreSegmentField(
        array.getUpdate(rowIndex + 1),
        Hdf5BottomSegment::genomeIndexOffset, end);
    row = ramaxReacquireSegmentRow(array, rowIndex, row);

    const hal_index_t nullIndex = NULL_INDEX;
    for (hal_size_t child = 0;
         child < _numChildren; ++child) {
        const size_t childOffset =
            Hdf5BottomSegment::firstChildOffset +
            child * (sizeof(hal_index_t) + sizeof(bool));
        ramaxStoreSegmentField(row, childOffset, nullIndex);
        const bool childReversed = false;
        ramaxStoreSegmentField(
            row, childOffset + sizeof(hal_index_t), childReversed);
    }
    ramaxStoreSegmentField(
        row, Hdf5BottomSegment::topIndexOffset, nullIndex);
}

void hal::RamaxHdf5SegmentWriter::setTopParentIndex(
    hal_index_t arrayIndex, hal_index_t parentIndex) {
    char *row = _genome->_topArray.getUpdate(
        static_cast<hsize_t>(arrayIndex));
    ramaxStoreSegmentField(
        row, Hdf5TopSegment::parentIndexOffset, parentIndex);
}

void hal::RamaxHdf5SegmentWriter::setTopBottomParseIndex(
    hal_index_t arrayIndex, hal_index_t bottomIndex) {
    char *row = _genome->_topArray.getUpdate(
        static_cast<hsize_t>(arrayIndex));
    ramaxStoreSegmentField(
        row, Hdf5TopSegment::bottomIndexOffset, bottomIndex);
}

void hal::RamaxHdf5SegmentWriter::setTopNextParalogyIndex(
    hal_index_t arrayIndex, hal_index_t nextIndex) {
    assert(nextIndex != arrayIndex);
    char *row = _genome->_topArray.getUpdate(
        static_cast<hsize_t>(arrayIndex));
    ramaxStoreSegmentField(
        row, Hdf5TopSegment::parIndexOffset, nextIndex);
}

void hal::RamaxHdf5SegmentWriter::setBottomTopParseIndex(
    hal_index_t arrayIndex, hal_index_t topIndex) {
    assert(arrayIndex >= 0);
    char *row = _genome->_bottomArray.getUpdate(
        static_cast<hsize_t>(arrayIndex));
    ramaxStoreSegmentField(
        row, Hdf5BottomSegment::topIndexOffset, topIndex);
}

void hal::RamaxHdf5SegmentWriter::setBottomChild(
    hal_index_t arrayIndex, hal_size_t child, hal_index_t childIndex,
    bool childReversed) {
    assert(arrayIndex >= 0);
    assert(child < _numChildren);
    char *row = _genome->_bottomArray.getUpdate(
        static_cast<hsize_t>(arrayIndex));
    const size_t childOffset =
        Hdf5BottomSegment::firstChildOffset +
        child * (sizeof(hal_index_t) + sizeof(bool));
    ramaxStoreSegmentField(row, childOffset, childIndex);
    ramaxStoreSegmentField(
        row, childOffset + sizeof(hal_index_t), childReversed);
}
]=])

set(ramax_bulk_dna_header [=[#ifndef _RAMAX_HDF5_BULK_DNA_H
#define _RAMAX_HDF5_BULK_DNA_H

#include "halDefs.h"

namespace hal {
class Genome;
class Hdf5Genome;

// Write ASCII DNA to the native packed HDF5 array without the per-base
// DnaIterator dispatch. The range is in genome coordinates. Ambiguous IUPAC
// input follows the historical cactus2hal policy: preserve case and map to N.
void ramaxWriteBulkDna(Genome &genome, hal_index_t start, const char *dna,
                       hal_size_t length);

// Bind once to a concrete HDF5 genome, then update complete segment rows
// without per-field virtual dispatch or repeated active-page discovery.
class RamaxHdf5SegmentWriter {
  public:
    explicit RamaxHdf5SegmentWriter(Genome &genome);

    void initializeTop(hal_index_t arrayIndex, hal_index_t start,
                       hal_size_t length, bool parentReversed);
    void initializeBottom(hal_index_t arrayIndex, hal_index_t start,
                          hal_size_t length);
    void setTopParentIndex(hal_index_t arrayIndex, hal_index_t parentIndex);
    void setTopBottomParseIndex(hal_index_t arrayIndex,
                                hal_index_t bottomIndex);
    void setTopNextParalogyIndex(hal_index_t arrayIndex,
                                 hal_index_t nextIndex);
    void setBottomTopParseIndex(hal_index_t arrayIndex,
                                hal_index_t topIndex);
    void setBottomChild(hal_index_t arrayIndex, hal_size_t child,
                        hal_index_t childIndex, bool childReversed);

  private:
    Hdf5Genome *_genome;
    hal_size_t _numChildren;
};
}

#endif
]=])
set(metadata_source
    "${RAMAX_HAL_SOURCE_DIR}/api/hdf5_impl/hdf5MetaData.cpp")
_ramax_read_pinned_hal_source(
    "${metadata_source}"
    "ab792be92884f001e169b9068af3e25abf801d419139fac4cb2ef77f9981f718"
    metadata_contents)
set(before_metadata_destructor [=[HDF5MetaData::~HDF5MetaData() {
    write();
}]=])
set(after_metadata_destructor [=[HDF5MetaData::~HDF5MetaData() {
    try {
        write();
    } catch (...) {
    }
}]=])
_ramax_replace_hal_block(metadata_contents before_metadata_destructor
    after_metadata_destructor "metadata non-throwing destruction")

set(array_source
    "${RAMAX_HAL_SOURCE_DIR}/api/hdf5_impl/hdf5ExternalArray.cpp")
_ramax_read_pinned_hal_source(
    "${array_source}"
    "25eef21d9223ada3815c01aac07e11c03a4b90d43cf338c3113305d547e7bacb"
    array_contents)
set(before_array_includes [=[#include <cassert>
#include <iostream>]=])
set(after_array_includes [=[#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>
#include <omp.h>
#include <zlib.h>

namespace {
// Compression never calls HDF5 from a worker. In-flight compressed bytes and
// the final padded chunk share this budget; at most eight zlib states are live.
constexpr size_t ramaxCompressionBufferBytes = 32U * 1024U * 1024U;
constexpr int ramaxCompressionMaxThreads = 8;
constexpr size_t ramaxCompressionMaxSlots = 256;

struct RamaxDeflater {
    z_stream stream{};
    bool initialized = false;

    RamaxDeflater() = default;
    RamaxDeflater(const RamaxDeflater&) = delete;
    RamaxDeflater& operator=(const RamaxDeflater&) = delete;
    ~RamaxDeflater() {
        if (initialized) {
            deflateEnd(&stream);
        }
    }
};

bool ramaxWriteDeflateChunks(H5::DataSet& dataset, const H5::DataType& type,
                            const char* input, hsize_t elements,
                            size_t elementBytes) {
#if !H5_VERSION_GE(1, 10, 2)
    return false;
#else
    const int requestedThreads =
        std::min(omp_get_max_threads(), ramaxCompressionMaxThreads);
    if (requestedThreads < 2 || omp_in_parallel() || elements == 0) {
        return false;
    }
    H5::DSetCreatPropList props = dataset.getCreatePlist();
    // Other filters, conversions and layouts retain HDF5's complete pipeline.
    if (props.getLayout() != H5D_CHUNKED || props.getNfilters() != 1 ||
        H5Tequal(type.getId(), dataset.getDataType().getId()) <= 0 ||
        H5Tdetect_class(type.getId(), H5T_VLEN) != 0 ||
        H5Tdetect_class(type.getId(), H5T_REFERENCE) != 0 ||
        (H5Tget_class(type.getId()) != H5T_STRING &&
         H5Tdetect_class(type.getId(), H5T_STRING) != 0) ||
        H5Tis_variable_str(type.getId()) != 0) {
        return false;
    }
    unsigned flags = 0;
    unsigned level = 0;
    size_t parameters = 1;
    const H5Z_filter_t filter = H5Pget_filter2(
        props.getId(), 0, &flags, &parameters, &level, 0, nullptr, nullptr);
    if (filter != H5Z_FILTER_DEFLATE || parameters != 1 || level > 9) {
        return false;
    }
    unsigned chunkOptions = 0;
    if (H5Pget_chunk_opts(props.getId(), &chunkOptions) < 0) {
        throw H5::DataSetIException("Hdf5ExternalArray::write",
                                   "Cannot read chunk filter options");
    }
    if (chunkOptions != 0) {
        return false;
    }
    hsize_t chunkElements = 0;
    if (props.getChunk(1, &chunkElements) != 1 || chunkElements == 0 ||
        elementBytes == 0 ||
        chunkElements > std::numeric_limits<uInt>::max() / elementBytes ||
        elements > std::numeric_limits<size_t>::max() / elementBytes) {
        return false;
    }
    const hsize_t chunks = elements / chunkElements +
                           (elements % chunkElements != 0);
    if (chunks < 2) {
        return false;
    }
    const size_t chunkBytes = static_cast<size_t>(chunkElements) * elementBytes;
    const uLong capacity = compressBound(static_cast<uLong>(chunkBytes));
    const bool partialLast = elements % chunkElements != 0;
    const size_t paddingBytes = partialLast ? chunkBytes : 0;
    if (capacity == 0 || capacity > std::numeric_limits<uInt>::max() ||
        paddingBytes >= ramaxCompressionBufferBytes) {
        return false;
    }
    const size_t slots = std::min<size_t>(
        std::min<hsize_t>(chunks, ramaxCompressionMaxSlots),
        (ramaxCompressionBufferBytes - paddingBytes) / capacity);
    if (slots < 2) {
        return false;
    }
    const int workers = std::min<int>(requestedThreads, static_cast<int>(slots));

    std::vector<unsigned char> lastChunk;
    if (partialLast) {
        H5D_fill_value_t fillStatus;
        if (H5Pfill_value_defined(props.getId(), &fillStatus) < 0) {
            throw H5::DataSetIException("Hdf5ExternalArray::write",
                                       "Cannot read chunk fill policy");
        }
        if (fillStatus == H5D_FILL_VALUE_UNDEFINED) {
            return false;
        }
        // Direct chunks must contain a full chunk, including the typed fill
        // outside the current extent. No uninitialized tail reaches the file.
        lastChunk.resize(chunkBytes);
        if (H5Pget_fill_value(props.getId(), type.getId(), lastChunk.data()) < 0) {
            throw H5::DataSetIException("Hdf5ExternalArray::write",
                                       "Cannot read chunk fill value");
        }
        for (size_t filled = elementBytes; filled < chunkBytes;) {
            const size_t count = std::min(filled, chunkBytes - filled);
            std::memcpy(lastChunk.data() + filled, lastChunk.data(), count);
            filled += count;
        }
        const size_t finalOffset =
            static_cast<size_t>(elements / chunkElements) * chunkBytes;
        std::memcpy(lastChunk.data(), input + finalOffset,
                    static_cast<size_t>(elements % chunkElements) * elementBytes);
    }

    std::vector<unsigned char> compressed(slots * capacity);
    std::unique_ptr<RamaxDeflater[]> deflaters(new RamaxDeflater[workers]);
    for (int worker = 0; worker < workers; ++worker) {
        const int result = deflateInit(&deflaters[worker].stream,
                                       static_cast<int>(level));
        if (result != Z_OK) {
            throw H5::DataSetIException("Hdf5ExternalArray::write",
                                       "Cannot initialize bounded zlib workers");
        }
        deflaters[worker].initialized = true;
    }
    std::array<int, ramaxCompressionMaxSlots> results{};
    std::array<uLong, ramaxCompressionMaxSlots> sizes{};
    for (hsize_t begin = 0; begin < chunks; begin += slots) {
        const int count = static_cast<int>(
            std::min<hsize_t>(slots, chunks - begin));
#pragma omp parallel for schedule(static) num_threads(workers)
        for (int slot = 0; slot < count; ++slot) {
            z_stream& stream = deflaters[omp_get_thread_num()].stream;
            int result = deflateReset(&stream);
            if (result == Z_OK) {
                const hsize_t chunk = begin + static_cast<hsize_t>(slot);
                const auto* bytes = partialLast && chunk == chunks - 1
                    ? lastChunk.data()
                    : reinterpret_cast<const unsigned char*>(input) +
                          static_cast<size_t>(chunk) * chunkBytes;
                stream.next_in = const_cast<Bytef*>(bytes);
                stream.avail_in = static_cast<uInt>(chunkBytes);
                stream.next_out = compressed.data() +
                                  static_cast<size_t>(slot) * capacity;
                stream.avail_out = static_cast<uInt>(capacity);
                result = deflate(&stream, Z_FINISH);
                if (result == Z_STREAM_END && stream.avail_in == 0) {
                    sizes[slot] = stream.total_out;
                    result = Z_OK;
                } else if (result == Z_STREAM_END || result == Z_OK) {
                    result = Z_BUF_ERROR;
                }
            }
            results[slot] = result;
        }
        // Check the complete batch before committing any of its chunks.
        for (int slot = 0; slot < count; ++slot) {
            if (results[slot] != Z_OK) {
                throw H5::DataSetIException("Hdf5ExternalArray::write",
                                           "Bounded zlib compression failed");
            }
        }
        for (int slot = 0; slot < count; ++slot) {
            const hsize_t offset =
                (begin + static_cast<hsize_t>(slot)) * chunkElements;
            if (H5Dwrite_chunk(dataset.getId(), H5P_DEFAULT, 0, &offset,
                              static_cast<size_t>(sizes[slot]),
                              compressed.data() +
                                  static_cast<size_t>(slot) * capacity) < 0) {
                throw H5::DataSetIException("Hdf5ExternalArray::write",
                                           "Cannot write compressed HAL chunk");
            }
        }
    }
    return true;
#endif
}
}  // namespace]=])
_ramax_replace_hal_block(array_contents before_array_includes after_array_includes
    "array initialization support")
set(before_array_allocation [=[    delete[] _buf;
    _buf = new char[_bufSize * _dataSize];]=])
set(after_array_allocation [=[    delete[] _buf;
    _buf = NULL;
    _dirty = false;
    _buf = new char[_bufSize * _dataSize];]=])
_ramax_replace_hal_block(array_contents before_array_allocation after_array_allocation
    "array allocation failure cleanup")
set(before_new_array_buffer [=[    // create the internal data buffer
    initBuf();

    // create the hdf5 array]=])
set(after_new_array_buffer [=[    // New arrays must not serialize previous heap contents in reserved
    // fields, padding, or end sentinels. Existing datasets are loaded by page().
    initBuf();
    std::memset(_buf, 0, _bufSize * _dataSize);

    // create the hdf5 array]=])
_ramax_replace_hal_block(array_contents before_new_array_buffer after_new_array_buffer
    "initialize newly created HAL arrays")
set(before_array_write [=[void Hdf5ExternalArray::write() {
    if (_dirty) {
        _dataSpace.selectHyperslab(H5S_SELECT_SET, &_bufSize, &_bufStart);
        _dataSet.write(_buf, _dataType, _chunkSpace, _dataSpace);
        _dirty = false;
    }
}]=])
set(after_array_write [=[void Hdf5ExternalArray::write() {
    if (_dirty) {
        if (_bufStart == 0 && _bufSize == _size &&
            ramaxWriteDeflateChunks(_dataSet, _dataType, _buf, _size, _dataSize)) {
            _dirty = false;
            return;
        }
        _dataSpace.selectHyperslab(H5S_SELECT_SET, &_bufSize, &_bufStart);
        _dataSet.write(_buf, _dataType, _chunkSpace, _dataSpace);
        _dirty = false;
    }
}]=])
_ramax_replace_hal_block(array_contents before_array_write after_array_write
    "bounded parallel compression with ordered HDF5 chunk writes")

set(sequence_source
    "${RAMAX_HAL_SOURCE_DIR}/api/hdf5_impl/hdf5Sequence.cpp")
_ramax_read_pinned_hal_source(
    "${sequence_source}"
    "a634a06d0bd7e9978dd2fcab934e5b095b10b9fac9b3a5f88bdb04e494f4632a"
    sequence_contents)
set(before_name_type [=[    StrType strType(PredType::NATIVE_CHAR, (maxNameLength + 1) * sizeof(char));]=])
set(after_name_type [=[    // NATIVE_CHAR is an integer datatype even inside StrType. Enlarging it
    // produces an invalid-width integer rejected by modern HDF5 readers.
    StrType strType(PredType::C_S1, (maxNameLength + 1) * sizeof(char));]=])
_ramax_replace_hal_block(sequence_contents before_name_type after_name_type
    "use a real fixed-width string datatype for sequence names")

file(MAKE_DIRECTORY "${RAMAX_HAL_PATCH_OUTPUT_DIR}")

# Every HDF5 implementation translation unit, every private header reachable
# from one of those units, and the outside alignment factory are copied into
# one source-pinned directory. Quoted includes therefore stay inside this
# overlay (not the source tree), so Hdf5Genome and HDF5DnaAccess have one class
# definition throughout the replacement set. Direct objects for the complete
# set prevent the corresponding original libHal.a members from being extracted.
_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5BottomSegment.cpp"
    "87aaa54f19c1c1df93512d0f1695cac05e445d209657d5f72c2e8b9f608574c5")
_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5SequenceIterator.cpp"
    "d244dec735079d87b7dae60e6e07148b930c8fce72422e5fd7618e48729bce02")
_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5TopSegment.cpp"
    "00c83deb5ded54dafa23cc406d157efafc521d365d400cdd27db91e472421e6e")
_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5UDCFuseDriver.cpp"
    "c1da920fb82d0587bdc0f2d0c0f0bf2d1262395d4106f753e6ff74c322be07c7")
_ramax_copy_pinned_hal_source(
    "api/impl/halAlignmentInstance.cpp"
    "c5474a4727ac197206672d1c574613cef2c0a0bc2112074f13ce9ac70eb50e9d")
_ramax_copy_pinned_hal_source(
    "api/impl/halCLParser.cpp"
    "69d602397ff7cfd50fc3ea65b67d7b63cb37c4d865be17e64baaf797deecb2e1")

_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5Alignment.h"
    "eb6cf4d5cf746234b51edf332ae55d08c6fe627c22d0aebc920dd8a37b1fecf7")
_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5Common.h"
    "82088aa1115eb40f2023ce037f69ac75a7876637189b74a6cf1726654cf277df")
_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5DnaArray.h"
    "5913b7f526c08d24570a6b008faf0b0457446d17d5d335382046ed0d29e2dd25")
_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5ExternalArray.h"
    "e65bcaf776bb6a7ae0b29d462b7cccdddc25992104ec0c636242d35a6a965d75")
_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5MetaData.h"
    "b43a3ec0877e0cfb908acc6c9bb3f341dd148e9955c60b782503092b70a50cc6")
_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5Sequence.h"
    "c7af6a4fa8a32b1825c740d29ed29a0d6f574501ad041bd9e710ba0e4e356adb")
_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5SequenceIterator.h"
    "a9dec286b9d7619df9e106106c87b8ad1bd48f8d2319122a49e38f80cc97872b")
_ramax_copy_pinned_hal_source(
    "api/hdf5_impl/hdf5UDCFuseDriver.h"
    "18e0f80191000a0598d8384b3b57a27c48908e82ea218edac7687fbbff5375f6")
file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/hdf5Alignment.cpp"
    "${alignment_contents}")
file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/hdf5Genome.cpp"
    "${genome_contents}")
file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/hdf5Genome.h"
    "${genome_header_contents}")
file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/hdf5BottomSegment.h"
    "${bottom_segment_header_contents}")
file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/hdf5TopSegment.h"
    "${top_segment_header_contents}")
file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/hdf5DnaDriver.cpp"
    "${dna_driver_contents}")
file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/hdf5DnaDriver.h"
    "${dna_driver_header_contents}")
file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/ramaxHdf5BulkDna.h"
    "${ramax_bulk_dna_header}")
file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/hdf5MetaData.cpp"
    "${metadata_contents}")
file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/hdf5ExternalArray.cpp"
    "${array_contents}")
file(WRITE "${RAMAX_HAL_PATCH_OUTPUT_DIR}/hdf5Sequence.cpp"
    "${sequence_contents}")
