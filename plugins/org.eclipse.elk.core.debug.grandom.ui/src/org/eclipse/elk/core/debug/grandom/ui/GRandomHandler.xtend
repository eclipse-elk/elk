/*******************************************************************************
 * Copyright (c) 2016, 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.grandom.ui

import java.util.Collections
import org.eclipse.core.commands.AbstractHandler
import org.eclipse.core.commands.ExecutionEvent
import org.eclipse.core.commands.ExecutionException
import org.eclipse.core.resources.IFile
import org.eclipse.core.runtime.IProgressMonitor
import org.eclipse.core.runtime.IStatus
import org.eclipse.core.runtime.Status
import org.eclipse.core.runtime.jobs.Job
import org.eclipse.elk.core.debug.grandom.gRandom.RandGraph
import org.eclipse.emf.common.util.URI
import org.eclipse.emf.ecore.resource.impl.ResourceSetImpl
import org.eclipse.jface.viewers.IStructuredSelection
import org.eclipse.ui.PlatformUI
import org.eclipse.ui.handlers.HandlerUtil
import org.eclipse.ui.statushandlers.StatusManager

class GRandomHandler extends AbstractHandler {

    override execute(ExecutionEvent event) throws ExecutionException {
        PlatformUI.getWorkbench().saveAllEditors(true);
        val selection = HandlerUtil.getCurrentSelection(event);
        if (selection instanceof IStructuredSelection) {
            val element = selection.getFirstElement();

            if (element instanceof IFile) {
                try {
                    // load the file into a resource  
                    val resourceSet = new ResourceSetImpl();
                    val r = resourceSet.createResource(URI.createFileURI("dummy.elkr"));
                    r.load(element.contents, Collections.emptyMap());
                    val project = element.project

                    if (!r.getContents().isEmpty() && r.getContents().get(0) instanceof RandGraph) {

                        val rand = r.contents.head as RandGraph

                        val job = new Job("Generate Random Graphs") {
                            override protected run(IProgressMonitor monitor) {
                                monitor.beginTask("Generate Random Graphs", 1);
                                try {
                                    new GRandomGraphMaker(rand).gen(project)
                                } catch (Exception e) {
                                    val status = new Status(IStatus.ERROR, "org.eclipse.elk.core.debug.grandom.ui", 0,
                                        "Random Graph Generation failed", e);
                                    StatusManager.getManager().handle(
                                        status,
                                        StatusManager.SHOW.bitwiseOr(StatusManager.LOG)
                                    );
                                } finally {
                                    monitor.done();
                                }
                                return Status.OK_STATUS;
                            }
                        };
                        job.setUser(true);
                        job.schedule();

                    }
                } catch (Exception e) {
//                    StatusManager.getManager().handle(
//                            new Status(IStatus.ERROR, GranaPlugin.PLUGIN_ID, MESSAGE_BATCH_FAILED, e));
                }
            }
        }

        return null;
    }

}
