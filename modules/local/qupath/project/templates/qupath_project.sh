#!/usr/bin/env bash

set -euo pipefail

echo '165e27a0731d58ba039e9d0d34a54acf896bb92e81922370df295cb85418ce53  QuPath-Linux.tar.xz' | sha256sum --check
tar -xJf QuPath-Linux.tar.xz QuPath/lib/app

cat > create_project.groovy <<'GROOVY'
import com.google.gson.JsonParser
import java.awt.image.BufferedImage
import java.nio.file.Files
import qupath.lib.images.ImageData
import qupath.lib.images.servers.ImageServerProvider
import qupath.lib.io.PathIO
import qupath.lib.projects.Projects

def directory = new File('qupath').getCanonicalFile()
def images = new File(directory, 'images')
def manifests = new File(directory, 'layers')
images.mkdirs()
manifests.mkdirs()
def project = Projects.createProject(new File(directory, 'project.qpproj'), BufferedImage.class)
def inputs = new File('slides').listFiles().sort { it.name }
def classes = new LinkedHashSet()
for (def input : inputs) {
    def sample = input.name
    def manifest = new File(input, 'layers.json')
    def expected = manifest.withReader { reader ->
        JsonParser.parseReader(reader).getAsJsonArray().sum { layer -> layer.get('objects').getAsLong() } ?: 0L
    }
    Files.copy(manifest.toPath(), new File(manifests, sample + '.json').toPath())

    def objectFiles = input.listFiles().findAll { it.name.endsWith('.geojson') }.sort { it.name }
    for (def slide : input.listFiles().findAll { it.name.endsWith('.ome.tif') }.sort { it.name == 'original.ome.tif' ? 0 : 1 }) {
        def target = new File(images, sample + '.' + slide.name)
        Files.copy(slide.toPath(), target.toPath())
        def support = ImageServerProvider.getPreferredUriImageSupport(
            BufferedImage.class, target.toURI().toString(), '--classname', 'BioFormatsServerBuilder'
        )
        def entry = project.addImage(support.getBuilders().get(0))
        entry.setImageName(slide.name == 'original.ome.tif' ? sample : sample + ' (normalized)')
        def data = entry.readImageData()
        long count = 0
        try {
            data.setImageType(ImageData.ImageType.BRIGHTFIELD_H_E)
            for (def file : objectFiles) {
                def objects = PathIO.readObjects(file.toPath())
                count += objects.size()
                objects.each { object ->
                    if (object.getPathClass() != null) {
                        classes.add(object.getPathClass())
                    }
                }
                data.getHierarchy().addObjects(objects)
            }
            if (count != expected || data.getHierarchy().nObjects() != expected) {
                throw new IllegalStateException('Object count mismatch for ' + entry.getImageName())
            }
            entry.saveImageData(data)
        } finally {
            data.getServer().close()
        }
    }
}
project.setPathClasses(new ArrayList(classes))
project.syncChanges()
GROOVY

java -Djava.awt.headless=true -Duser.home="\$PWD" \
    -XX:ActiveProcessorCount=${task.cpus} -Xmx${(task.memory.toMega() * 0.75).toInteger()}m \
    -cp 'QuPath/lib/app/*' qupath.QuPath script create_project.groovy

test -s qupath/project.qpproj
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    qupath: "0.7.0"
END_VERSIONS
