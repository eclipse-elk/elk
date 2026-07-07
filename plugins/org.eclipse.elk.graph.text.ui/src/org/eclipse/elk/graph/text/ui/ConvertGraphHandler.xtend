/*******************************************************************************
 * Copyright (c) 2017, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *     Miro Spönemann - converted to Xtend
 *******************************************************************************/
package org.eclipse.elk.graph.text.ui

import com.google.inject.Inject
import java.io.IOException
import java.util.Collections
import java.util.List
import org.eclipse.core.commands.AbstractHandler
import org.eclipse.core.commands.ExecutionEvent
import org.eclipse.core.commands.ExecutionException
import org.eclipse.core.resources.IFile
import org.eclipse.core.runtime.IPath
import org.eclipse.core.runtime.IProgressMonitor
import org.eclipse.core.runtime.IStatus
import org.eclipse.core.runtime.Path
import org.eclipse.core.runtime.Status
import org.eclipse.core.runtime.jobs.Job
import org.eclipse.elk.core.data.LayoutMetaDataService
import org.eclipse.elk.core.data.LayoutOptionData
import org.eclipse.elk.core.options.CoreOptions
import org.eclipse.elk.core.service.util.ProgressMonitorAdapter
import org.eclipse.elk.core.util.IElkProgressMonitor
import org.eclipse.elk.core.validation.LayoutOptionValidator
import org.eclipse.elk.graph.EMapPropertyHolder
import org.eclipse.elk.graph.ElkEdge
import org.eclipse.elk.graph.ElkGraphElement
import org.eclipse.elk.graph.ElkLabel
import org.eclipse.elk.graph.ElkNode
import org.eclipse.elk.graph.ElkPort
import org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl
import org.eclipse.elk.graph.properties.IProperty
import org.eclipse.elk.graph.properties.IPropertyValueProxy
import org.eclipse.elk.graph.util.GraphIdentifierGenerator
import org.eclipse.emf.common.util.URI
import org.eclipse.emf.ecore.EObject
import org.eclipse.emf.ecore.resource.impl.ResourceSetImpl
import org.eclipse.emf.ecore.util.EcoreUtil
import org.eclipse.jface.viewers.IStructuredSelection
import org.eclipse.ui.handlers.HandlerUtil
import org.eclipse.ui.statushandlers.StatusManager

import static extension org.eclipse.elk.graph.text.ElkGraphTextUtil.*
import static extension org.eclipse.xtext.EcoreUtil2.*

/** 
 * A command handler that can be used to convert a graph from the textual format to the XMI format
 * and back. The default implementation only copies input models and makes certain modifications if
 * the target format is the textual format.
 * 
 * <p>Note that plug-ins containing subclasses of this handler must define their own convert graphs
 * command. Only then can they declare menu contributions that reference their new command, which
 * in turn defines their new subclass as command handler.</p>
 */
class ConvertGraphHandler extends AbstractHandler {
    
    /** Parameter id for the target file extension. */
    static final String PARAM_TARGET_EXT = 'org.eclipse.elk.graph.text.convert.targetExtension'
    /** Textual graph file extension. */
    static final String EXT_ELK_TEXT = 'elkt'
    /** XMI graph file extension. */
    static final String EXT_ELK_XMI = 'elkg'
    
    /** The target file extension. */
    String targetExtension
    
    @Inject
    LayoutOptionValidator layoutOptionValidator

    override final Object execute(ExecutionEvent event) throws ExecutionException {
        targetExtension = event.getParameter(PARAM_TARGET_EXT)
        if (!targetExtension.isValidTargetExtension) {
            val status = new Status(IStatus.ERROR,
                ElkGraphUiModule.PLUGIN_ID, '''Invalid target extension: «targetExtension»''')
            StatusManager.manager.handle(status, StatusManager.SHOW + StatusManager.LOG)
            return null
        }
        val selection = HandlerUtil.getCurrentSelection(event)
        if (selection instanceof IStructuredSelection) {
            val files = selection.toList.filter(IFile)
            
            // Put the conversion into a job with a progress bar,
            // because the conversion may take a while :-/ 
            val Job job = new Job("Convert models") {
                override protected run(IProgressMonitor monitor) {
                    val elkMonitor = new ProgressMonitorAdapter(monitor)
                    elkMonitor.begin("Convert models", files.size)
                    files.forEach[ file |
                        if (!elkMonitor.isCanceled)
                            file.convert(elkMonitor.subTask(1))
                    ]
                    elkMonitor.done
                    return Status.OK_STATUS
                }
            }
            job.user = true
            job.schedule()
        }
        return null
    }

    private def isValidTargetExtension(String targetExtension) {
        targetExtension == EXT_ELK_TEXT || targetExtension == EXT_ELK_XMI
    }

    /** 
     * Process a source file.
     */
    def void convert(IFile file, IElkProgressMonitor monitor) {
        monitor.begin('''Convert «file.name»''', 3)
        try {
            val sourceModels = readModel(URI.createPlatformResourceURI(file.fullPath.toString, false))
            monitor.worked(1)
            val targetModels = (sourceModels as Iterable<EObject>).map[transform].toList
            monitor.worked(1)
            saveModel(targetModels, createTarget(file))
        } catch (Exception exception) {
            val status = new Status(IStatus.ERROR, ElkGraphUiModule.PLUGIN_ID,
                "Error while converting the selected graph.", exception)
            StatusManager.manager.handle(status, StatusManager.SHOW + StatusManager.LOG)
        } finally {
            monitor.done
        }
    }

    /** 
     * Transform the graph before it is written to the new file format.
     */
    protected def EObject transform(EObject model) {
        val copy = EcoreUtil.copy(model)
        if (targetExtension.equals(EXT_ELK_TEXT) && copy instanceof ElkNode) {
            val graph = copy as ElkNode
            // We want to convert to the textual format, so write missing identifiers into the graph
            GraphIdentifierGenerator.forGraph(graph)
                .assertValid
                .assertExists
                .assertUnique
                .execute()
            // Uncomment the following lines to reduce the graph to what's really relevant
//             removeInvalidProperties(graph)
//             removeUnnecessaryLayouts(graph)
        }
        return copy
    }

    /** 
     * Create a target file URI from a source file.
     */
    private def createTarget(IFile file) {
        val basePath = file.fullPath
        val name = basePath.removeFileExtension.lastSegment
        var IPath targetPath = new Path('''«name».«targetExtension»''')
        var i = 0
        val parent = file.parent
        while (parent.exists(targetPath)) {
            targetPath = new Path('''«name»«i++».«targetExtension»''')
        }
        targetPath = parent.fullPath.append(targetPath)
        return URI.createPlatformResourceURI(targetPath.toString, false)
    }

    /** 
     * Read the model contents from the given URI.
     */
    private def readModel(URI uri) {
        val resourceSet = new ResourceSetImpl
        val resource = resourceSet.getResource(uri, true)
        return resource.contents
    }

    /** 
     * Save a collection of models to the given URI.
     */
    private def saveModel(List<EObject> models, URI uri) throws IOException {
        val resourceSet = new ResourceSetImpl
        val resource = resourceSet.createResource(uri)
        resource.contents.addAll(models)
        resource.save(Collections.EMPTY_MAP)
    }
    
    protected def removeInvalidProperties(ElkNode graph) {
        val toRemove = newArrayList
        graph.eAllContents.filter(EMapPropertyHolder).forEach[ element |
            element.properties.forEach[ entry |
                if (!isValid(entry.key, entry.value, element))
                    toRemove += entry as ElkPropertyToValueMapEntryImpl
            ]
        ]
        toRemove.forEach[EcoreUtil.remove(it)]
    }
    
    static val EXCLUDED_OPTIONS = #[
        CoreOptions.POSITION, CoreOptions.BEND_POINTS, CoreOptions.JUNCTION_POINTS
    ]
    
    private def isValid(IProperty<?> property, Object value, EMapPropertyHolder container) {
        // Check whether the property is a registered layout option
        val option =
            if (property instanceof LayoutOptionData)
                property
            else
                LayoutMetaDataService.instance.getOptionData(property.id)
        if (option === null)
            return false
        if (EXCLUDED_OPTIONS.contains(option))
            return false
        
        // Check whether the value is valid
        val parsedValue =
            if (value instanceof String)
                option.parseValue(value)
            else if (value instanceof IPropertyValueProxy)
                value.resolveValue(option)
            else
                value
        if (parsedValue === null)
            return false
        
        if (container instanceof ElkGraphElement) {
            // Check whether the option target is valid
            switch container {
                ElkNode:
                    if (!(option.targets.contains(LayoutOptionData.Target.NODES)
                            || !container.children.empty && option.targets.contains(LayoutOptionData.Target.PARENTS)))
                        return false
                ElkEdge:
                    if (!option.targets.contains(LayoutOptionData.Target.EDGES))
                        return false
                ElkPort:
                    if (!option.targets.contains(LayoutOptionData.Target.PORTS))
                        return false
                ElkLabel:
                    if (!option.targets.contains(LayoutOptionData.Target.LABELS))
                        return false
            }
            
            // Check whether the option is supported by the layout algorithm
            if (container instanceof ElkNode && option.targets.contains(LayoutOptionData.Target.NODES)) {
                if (!isSupportedByAlgorithm(option, (container as ElkNode).parent))
                    return false
            } else {
                if (!isSupportedByAlgorithm(option, container))
                    return false
            }
            
            // Check for any issues reported by the layout option validator
            val issues = layoutOptionValidator.checkProperty(option, parsedValue, container)
            if (!issues.empty)
                return false
            
            return true
        } else {
            return false
        }
    }
    
    private def isSupportedByAlgorithm(LayoutOptionData option, ElkGraphElement element) {
        if (option == CoreOptions.ALGORITHM)
            return true
        var algorithm = element.algorithm
        if (algorithm !== null && algorithm.knowsOption(option))
            return true
        // Maybe the option is handled by one of the parent nodes
        var parent = element.getContainerOfType(ElkNode)
        while (parent !== null) {
            algorithm = parent.algorithm
            if (algorithm !== null && algorithm.knowsOption(option))
                return true
            parent = parent.parent
        }
        return false
    }
    
    protected def removeUnnecessaryLayouts(ElkNode graph) {
        val toRemove = newArrayList
        graph.eAllContents.filter(ElkEdge).forEach[ edge |
            toRemove += edge.sections
        ]
        toRemove.forEach[EcoreUtil.remove(it)]
    }
    
}
